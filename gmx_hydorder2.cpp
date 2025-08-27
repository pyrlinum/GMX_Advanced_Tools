// This file is a fork of GROMACS ver 2025 GMX HYDORDER TOOL
// Changes with respet to orignal code:
// This version computes frequency distribution of OTO parameter in either original or normalized form
// over the whole box or a volume bewteen surfaces defined by lower and upper radius from 'Protein' group
// 
// This is an unoptimized serial version 



#include "gmxpre.h"

#include <cmath>
#include <cstdio>
#include <cstring>

#include <filesystem>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/rgb.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/binsearch.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/gmxana/powerspect.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

// Define group names:
enum class PbcType : int;
struct gmx_output_env_t;

#define Q4DEF -1000.0

//-----------------------------------------------------------------------------------------------------------------------------
// Utlity inline functions:
// Alternative definitions of Orientational Tetrahedral Order Parameter:
// (form = 0) - original Chau-Hardwick(1998) parameter: 0 - tetrahedral order, 0.25 - random distribution
// (form = 1) - rescaled Errington-Debenedetti(2001) parameter: 1 - tetrahedral order, 0 - ideal gas. Can occasonally give produce values
inline real gmx_hyd2_oto_func(real x,int q4form) { return (q4form==0)  ? 3./32.*x :  1.-3./8.*x ;}
//============================================================================================================================
// Auxhiliary Declarations:
//============================================================================================================================
// Define groups and indices:
static void rd_groups_by_index(gmx::ArrayRef<const IndexGroup> indexGroups, int ngrps, int grpIdx[], char* gnames[], int isize[], int* index[]);
// Calculate OTO parameter per atom and save to qarr (if out of shell - -1)
static void calc_tetra_order_vals(t_topology top, PbcType  pbcType, int  natoms, matrix box,rvec x[],int ng,int* isize,int** index,int sflag,real rlo,real rhi,int q4form,real* qarr);

// Wrappers:
// Calculate and print out OTO values over trajectory:
static void calc_tetra_order_values(const char* fnTPS,const char* fnTRX,int ng,int* isize,int** index,int tblock,int* nframes,int sflag,real rlo,real rhi,int q4form,gmx_output_env_t* oenv);
// Accumulate OTO hystogram over trajectory:
static void calc_tetra_order_histogram(const char* fnTPS, const char* fnTRX,int ng,int* isize,int** index,int tblock,int* nframes,int sflag,real rlo,real rhi,int q4form,int hist_bins,int* sg_hist,real* sg_bin_centers, gmx_output_env_t* oenv);

// Helper functions:
// calculate 1D histrogram:
static void calc_hist( int maxidx, real* qarr, int hist_bins,real hist_bin_width, int* sg_hist);
// Output of OTO values per frame:
static void write_q4vals(int frame, int natm_write, int*  index, real* qarr);
// Output of hystogram:
static void writehist(int nframes,real* sg_bin_centers, int* sg_hist, int hist_bins, gmx::ArrayRef<const std::string> fnms);

//==========================================================================================================================
// Main program:
//==========================================================================================================================
int gmx_hydorder2(int argc, char* argv[])
{
    static const char* desc[] = {
        "[THISMODULE] is a fork of GROMACS HYDORDER module to compute the ",
        "tetrahedrality order parameter around a given atom type (defalut: O_WAT).",
        "Unlike the original, this module computs a histogram of the parameter",
        "time-averaged over selected tpart of the trajectory",
        "In this version only angle order parameter is calculated.",
        "Two versions of order formula are implemented and selected by -q4form argument:",
        "(0) - the original version by Chau and A.J. Hardwick (1998)",
        "(1) - the rescaled version by Errington and Debenedetti (2001)",
    };

    static const real EPS = 1.0e-6;
    static int   nsttblock = 1;
    static int   i, frames = 0;
    static int   q4bins = 0; 
    /* default: 0 - do not compute OTO parameter
    if 1 - compute individual OTO values per atom for each frame in separate file
    if >1 - compute a sigle histogram */
    static int   q4form = 1;
    /*   0 - Chau & Hardwick 1998 - Gromacs default: tatrahedral water q4=0, random bonds q4=0.25
    1 - Errington & Debenedetti 2001 - average value of q varies from 0 for an ideal gas to 1 for a regular tetrahedron */

    // define default group names:
    static const char* water_gname = "O_&_WAT";
    static const char* solute_gname = "Protein";
    static const char* ancor_gname = "O*_N*_&_Protein";
    static const char* exclude_gname = "Exclude";
    const char* grpDescr[] = {water_gname, solute_gname, ancor_gname, exclude_gname};

    // default output file names:
    static const char* ofile_hist = "oto_histogram.txt";

    // Cut-off radii used to split target atoms (water oxygens) by proximity to another group (protein)
    static real rcut_lower = -1.0;    // lower bound (<0.0 - no cutoff)
    static real rcut_upper = -1.0;    // upper bound (<0.0 - no cutoff )
    bool shell_flag = 0;             // whether to compute OTO over hydration shell only or the whole water system
    int hist_nbins;                  // number of bins to store histogram
    real wbin;                       // width of one histogram bin
    int *q4hist_values = nullptr;
    real *q4bin_centers = nullptr;

    // Groups to study: 1 - water, 2 - hydration shell 3 - additional atoms from solute to use as vertices in OTO calculation
    int grpNum = 4;
    int grpIdx[] = {-1,-1,-1,-1};
    char** grpname = nullptr;
    int**  index   = nullptr;
    int*   isize   = nullptr;


    // Numeric command line arguments:
    t_pargs pa[] = {
        { "-rlo", FALSE, etREAL, { &rcut_lower }, "Lower cut-off radius for hydration shell" },
        { "-rhi", FALSE, etREAL, { &rcut_upper }, "Upper cut-off radius for hydration shell" },
        { "-tblock", FALSE, etINT, { &nsttblock }, "Number of frames in one time-block average" },
        { "-q4bins", FALSE, etINT, { &q4bins}, "Number of bins in histogram of tetrahedral order parameter" },
        { "-q4form", FALSE, etINT, { &q4form}, "Form of tetrahedral order parameter: 0 - original Chau-Hardwick 1 - rescaled Errington-Debenedetti"        }
    };

    // Files for input and output:
    t_filenm fnm[] = {
        { efTRX, "-f", nullptr, ffREAD },              /* trajectory file              */
        { efNDX, "-n", nullptr, ffREAD },              /* index file           */
        { efTPR, "-s", nullptr, ffREAD },              /* topology file                */
        { efOUT, "-or", "raw", ffOPTWRMULT },          /* histogram output file         */
    };
#define NFILE asize(fnm)

    // Parsing arguments:-------------------------------------------------------------------------------------------------
    /*Filenames*/
    const char *      ndxfnm, *tpsfnm, *trxfnm;
    gmx_output_env_t* oenv;

    if (!parse_common_args(
                &argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    // inputs:
    ndxfnm = ftp2fn(efNDX, NFILE, fnm);
    tpsfnm = ftp2fn(efTPR, NFILE, fnm);
    trxfnm = ftp2fn(efTRX, NFILE, fnm);


    // Output file:
    gmx::ArrayRef<const std::string> raw = opt2fns("-or", NFILE, fnm);
    if ( (q4bins > 1 ) && (raw.size() != 1))
    {
        gmx_fatal(FARGS, "No or not correct number (1) of output-files: %td \n", raw.ssize());
    }

    // Number of the histogram bins:
    if (q4bins <= 0)
    {
        gmx_fatal(FARGS, "No or incorrect number of bins is defined: %d \n", q4bins);
    }

    // Formula for Orientational tetrahedral order parameter (q4):
    if ( q4form == 0 ) {
        fprintf(stdout,"Using the original Chau-Hardwick 1997 parameter\n");
    } else if ( q4form == 1 ) {
        fprintf(stdout,"Using the rescaled Errington-Debenedetti 2001 parameter\n");
    } else {
        gmx_fatal(FARGS, "Unexpected value (%d) for q4 functinal form: either 0 or 1 is implemented\n", q4form);
    }

    // Cut-offs sanity check:
    if ( ( (rcut_lower>=0.0) && (rcut_upper>0.0) ) && (rcut_lower > rcut_upper) )
    {
        gmx_fatal(FARGS, "Incorrect cutoffs for hydration shell are defined: rlo=%0.3f nm rhi=%0.3f nm \n",rcut_lower,rcut_upper );
    } else if ( (rcut_upper>0.0) || ((rcut_lower>0.0) && (rcut_upper==-1.0) )) {
        shell_flag = 1;
    }

    // Analyze groups:-------------------------------------------------------------------------------------------------
    std::vector<IndexGroup> ndxGrps = init_index(ndxfnm);
    fprintf(stdout,"===========================================================\n");
    fprintf(stdout,"Analysing Index file for predefined group names:\n");
    for(int ig=0;ig<grpNum;ig++) {
        grpIdx[ig] = find_group(grpDescr[ig],ndxGrps);
        if (grpIdx[ig]>=0) {
            printf("Group %s found as group index %d\n",grpDescr[ig],grpIdx[ig]);
        }
    }
    // Whether we are doing bulk orhydration shell claculation:
    if (grpIdx[0]==-1) {
         gmx_fatal(FARGS, "ERROR: Solvent molecules are not detected or group name does not match expected group name %s  \n",water_gname );
    } else if ( shell_flag && (grpIdx[0]==-1) ) {
        shell_flag = 0 ;
        rcut_lower = -1.0;
        rcut_upper = -1.0;
        printf("WARNING: Protein molecule not found - ignoring cut-off radii\n");
    }
    // allocate groups
    snew(isize,grpNum);
    snew(grpname,grpNum);
    snew(index,grpNum);
    rd_groups_by_index(ndxGrps,grpNum,grpIdx,grpname,isize,index);

    // Conclusion before calculation:
    fprintf(stdout,"===========================================================\n");
    if ( shell_flag == 0 ) {
        fprintf(stdout,"Bulk calculation: accumulationg order parameter for all water molecules\n" );
    } else {
        fprintf(stdout,"Shell calculation: accumulationg order parameter for the molecules between rlo=%0.3f nm rhi=%0.3f nm \n",rcut_lower,rcut_upper );
    }
    if ( isize[2]>0 ) {
        fprintf(stdout,"Group %s was detected - using additional atoms as vertices for OTO calculation\n",ancor_gname );
    }

    if ( isize[3]>0 ) {
        fprintf(stdout,"Group %s was detected - excluding from OTO calculation on each step those water molecules, which have at least one of nearest neighbours from this group\n",exclude_gname );
    }

   
    if ( q4bins > 1 ) {
        /* Calculate q4 histogram and print it out as raw output:*/
        // add first and last bins to collect out of [0:1] bounds values
        fprintf(stdout,"Calculating distribution of OTO values per each of the selected atoms\n");

        hist_nbins = q4bins+2;
        snew(q4hist_values,hist_nbins);
        snew(q4bin_centers,hist_nbins);
        
        // calculating global OTO distribution:
        calc_tetra_order_histogram(tpsfnm, trxfnm, grpNum, isize, index, nsttblock, &frames, shell_flag, rcut_lower, rcut_upper, q4form, hist_nbins, q4hist_values, q4bin_centers, oenv);

        writehist(frames,q4bin_centers, q4hist_values, hist_nbins, raw);

        // Free allocated memory:
        sfree(q4hist_values);
        sfree(q4bin_centers);
        for(int g=0;g<grpNum;g++) {
            if (isize[g]>0) { sfree(index[g]); }
            sfree(grpname[g]);
        }

    } else if ( q4bins == 1 ) {
        /* Calculate individual q4 values and print it out as raw output:*/
        fprintf(stdout,"Calculating individual OTO values per each of the selected atoms\n");
        calc_tetra_order_values(tpsfnm, trxfnm, grpNum, isize, index, nsttblock, &frames, shell_flag, rcut_lower, rcut_upper, q4form, oenv);           
    } 


    sfree(isize);
    sfree(index);
    sfree(grpname);

    return 0;
}

//=============================================================================================================================
// Auxhiliary functions:
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Main compute - calculate OTO in a frame and save to array qarr:
static void calc_tetra_order_vals(t_topology top,
                                  PbcType    pbcType,
                                  int        natoms,
                                  matrix     box,
                                  rvec       x[],
                                  int        ng,
                                  int*       isize,
                                  int**      index,
                                  int        sflag,
                                  real       rlo,
                                  real       rhi,
                                  int        q4form,
                                  real*      qarr)
{
    const int   maxidx = isize[0]; // number of water oxygens in total system
    const int   nei_grp[2] = { 0, 2 }; // predefined groups to scan for nearest neighbours
    int         ix, jx, i, j, k, gndx, ns, nxcl, *nn[4], *nng[4];
    int         slindex_h;
    bool        *smask, in_shell_flag;
    rvec        dx, rj, rk, urk, urj;
    real        cost, cost2, *sgmol, r2, rl2, rh2, box2, *r_nn[4];
    t_pbc       pbc;
    real        onethird = 1.0 / 3.0;
    gmx_rmpbc_t gpbc;

    box2 = box[XX][XX] * box[XX][XX];

    // allocate memory:
    for (i = 0; (i < 4); i++)
    {
        snew(r_nn[i], natoms);
        snew(nn[i], natoms);
        snew(nng[i], natoms);

        for (j = 0; (j < natoms); j++)
        {
            r_nn[i][j] = box2;
        }
    }
    snew(sgmol, maxidx); // array of tetrahedral order parameter
    snew(smask, maxidx); // mask array if include an atom into OTO calculation as center atom

    // Must init pbc every step because of pressure coupling
    set_pbc(&pbc, pbcType, box);
    gpbc = gmx_rmpbc_init(&top.idef, pbcType, natoms);
    gmx_rmpbc_apply(gpbc, natoms, box, x);

    rl2 = 0.0 ; 
    rh2 = box2 ;
    // Step 0. if shell calculation - select oxygen atoms within shell of protein:
    if (sflag) {
        // sanity check on rlo/rhi:
        rl2 =  (rlo>0) ? rlo*rlo : 0.0 ; 
        rh2 =  (rhi>0) ? rhi*rhi : box2 ;

        //TODO: Replace with 2 independent loops for each condition, then multiply matrixes.
        for (i = 0; (i < maxidx); i++)
        {
            ix = index[0][i];
            // First run - include all closer then rhi
            if ( rhi < 0 ) { // no upper limit:
                in_shell_flag = true;
            } else {
                in_shell_flag = false;
                for (j = 0; (j < isize[1]); j++)
                {
                    jx = index[1][j];
                    pbc_dx(&pbc, x[ix], x[jx], dx);
                    r2 = iprod(dx, dx);
                    in_shell_flag |= (r2<=rh2);
                }
            }

            // Second run - exclude all closer then rlo
            if ( rlo > 0 ) {
                for (j = 0; (j < isize[1]); j++)
                {
                    jx = index[1][j];
                    pbc_dx(&pbc, x[ix], x[jx], dx);
                    r2 = iprod(dx, dx);
                    in_shell_flag &= (r2>rl2);
                }
            }

            smask[i] = in_shell_flag;
        }
    } else { // bulk calculation - include all water oxygens
        (void)memset(smask, true, maxidx);
    }

    nxcl=0;
    // Main OTO Calculation:
    for (i = 0; (i < maxidx); i++)
        if (smask[i]) { // loop over oxygens inside shell (or all O group if no shell required):
            ix = index[0][i];

            // Step 1. Loop over neighbours: (0 - all water oxygens, 2 - additional ancor atoms if provided):
            for(int ig=0; ig<2; ig++) {
                gndx = nei_grp[ig];
                for (j = 0; (j < isize[gndx]); j++)
                {

                    if ((i == j) && (gndx == 0) ) // skip if the same atom of the same group
                    {
                        continue;
                    }

                    jx = index[gndx][j];

                    pbc_dx(&pbc, x[ix], x[jx], dx);
                    r2 = iprod(dx, dx);

                    // determine the nearest neighbours
                    if (r2 < r_nn[0][i])
                    {
                        r_nn[3][i] = r_nn[2][i];
                        nn[3][i]   = nn[2][i];
                        nng[3][i]  = nng[2][i];

                        r_nn[2][i] = r_nn[1][i];
                        nn[2][i]   = nn[1][i];
                        nng[2][i]  = nng[1][i];

                        r_nn[1][i] = r_nn[0][i];
                        nn[1][i]   = nn[0][i];
                        nng[1][i]  = nng[0][i];

                        r_nn[0][i] = r2;
                        nn[0][i]   = j;
                        nng[0][i]  = gndx;
                    }
                    else if (r2 < r_nn[1][i])
                    {
                        r_nn[3][i] = r_nn[2][i];
                        nn[3][i]   = nn[2][i];
                        nng[3][i]  = nng[2][i];

                        r_nn[2][i] = r_nn[1][i];
                        nn[2][i]   = nn[1][i];
                        nng[2][i]  = nng[1][i];

                        r_nn[1][i] = r2;
                        nn[1][i]   = j;
                        nng[1][i]  = gndx;
                    }
                    else if (r2 < r_nn[2][i])
                    {
                        r_nn[3][i] = r_nn[2][i];
                        nn[3][i]   = nn[2][i];
                        nng[3][i]  = nng[2][i];

                        r_nn[2][i] = r2;
                        nn[2][i]   = j;
                        nng[2][i]  = gndx;
                    }
                    else if (r2 < r_nn[3][i])
                    {
                        r_nn[3][i] = r2;
                        nn[3][i]   = j;
                        nng[3][i]  = gndx;
                    }
                }
            }

            // Step 1B. if excluded group is present - check if one of the nearest neighbours is from this group:
            gndx=3;
            j=0;
            in_shell_flag = smask[i];
            while ((j < isize[gndx]) && (in_shell_flag))
            {
                if (in_shell_flag) {
                    jx = index[gndx][j];

                    pbc_dx(&pbc, x[ix], x[jx], dx);
                    r2 = iprod(dx, dx);

                    // determine if excluded atom is closer than other nearest neighbours
                    if (r2 < r_nn[3][i]+0.000001)
                    {
                        in_shell_flag=false;
                        smask[i] = false;
                        nxcl++;
                    }
                }
                j++;
            }

            if (smask[i])  {
                // Step 2. Calculate angular part tetrahedrality order parameter per atom
                sgmol[i] = 0.0;
                for (j = 0; (j < 3); j++)
                {
                    for (k = j + 1; (k < 4); k++)
                    {
                        pbc_dx(&pbc, x[ix], x[index[nng[k][i]][nn[k][i]]], rk); // take in account the group of neighbour atom
                        pbc_dx(&pbc, x[ix], x[index[nng[j][i]][nn[j][i]]], rj); // take in account the group of neighbour atom

                        unitv(rk, urk);
                        unitv(rj, urj);

                        cost  = iprod(urk, urj) + onethird;
                        cost2 = cost * cost;

                        sgmol[i] += cost2;
                    }
                }
                sgmol[i] = gmx_hyd2_oto_func(sgmol[i],q4form);
            }

            // Save tuples of neqighbours and q4 values to output arrays:
            if (smask[i])  {
                qarr[i] = sgmol[i];
            }
        } // end if (smask[i])

    sfree(smask);
    sfree(sgmol);
    for (i = 0; (i < 4); i++)
    {
        sfree(r_nn[i]);
        sfree(nn[i]);
        sfree(nng[i]);
    }
}
//-----------------------------------------------------------------------------------------------------------------------------
// Wrappers:
// Calculate and print out OTO values over trajectory:
static void calc_tetra_order_values(const char*       fnTPS,
                                       const char*       fnTRX,
                                       int               ng,
                                       int*              isize,
                                       int**             index,
                                       int               tblock,
                                       int*              nframes,
                                       int               sflag,
                                       real              rlo,
                                       real              rhi,
                                       int               q4form,
                                       gmx_output_env_t* oenv)
{
    t_topology   top;
    PbcType      pbcType;
    t_trxstatus* status;
    int          natoms;
    real         t,binw;
    rvec *       xtop, *x;
    matrix       box;
    real*        qarr = nullptr; 
    int          i, framenr;
    const real onehalf = 1.0 / 2.0;


    // Check topology and trajectory files for consistency with index file:
    read_tps_conf(fnTPS, &top, &pbcType, &xtop, nullptr, box, FALSE);
    // check if atom numbers in Topology and Trajectory files match:
    natoms = read_first_x(oenv, &status, fnTRX, &t, &x, box);
    if (natoms > top.atoms.nr) { gmx_fatal(FARGS, "Topology (%d atoms) does not match trajectory (%d atoms)", top.atoms.nr, natoms); }
    // check if all atoms in selected index groups matching trajectory atoms:
    for(int ig=0; ig<ng; ig++) { if (isize[ig]>0) { check_index(nullptr, isize[ig], index[ig], nullptr, natoms); }}

    snew(qarr,isize[0]);
    for(int iq=0;iq<isize[0];iq++) qarr[iq]=Q4DEF;
    *nframes = 0;
    framenr  = 0;

    /* Loop over frames*/
    //TODO: add time average - accumulate x[] from NT read_next_x() and pass x_avg and box_avg to calc_tetra_order_hist()
    fprintf(stdout, "Processing frames:\n");
    do
    {
        // sample histogram every tblock steps (use the first frame without averaging for consistency)
        // Can calculate running average of the atom coordinates here before main calculator
        if (framenr % tblock == 0)
        {
            
            for(int iq=0;iq<isize[0];iq++) qarr[iq]=Q4DEF;
            
            calc_tetra_order_vals(top, pbcType, natoms, box, x, ng, isize, index, sflag, rlo, rhi, q4form, qarr);
            // save to output file:
            write_q4vals(framenr, isize[0], index[0], qarr);

            // update collected frame count:
            (*nframes)++;

        }
        // update total frame count:
        framenr++;

    } while (read_next_x(oenv, status, &t, x, box));
    close_trx(status);
    sfree(qarr);
}

// Accumulate OTO hystogram over trajectory:
static void calc_tetra_order_histogram(const char*       fnTPS,
                                       const char*       fnTRX,
                                       int               ng,
                                       int*              isize,
                                       int**             index,
                                       int               tblock,
                                       int*              nframes,
                                       int               sflag,
                                       real              rlo,
                                       real              rhi,
                                       int               q4form,
                                       int               hist_bins,
                                       int*              sg_hist,
                                       real*             sg_bin_centers,
                                       gmx_output_env_t* oenv)
{
    t_topology   top;
    PbcType      pbcType;
    t_trxstatus* status;
    int          natoms;
    real         t,binw;
    rvec *       xtop, *x;
    matrix       box;
    int*         sg_hist_tblk = nullptr;
    real*        qarr = nullptr; 
    int          i, framenr;
    const real onehalf = 1.0 / 2.0;


    // Check topology and trajectory files for consistency with index file:
    read_tps_conf(fnTPS, &top, &pbcType, &xtop, nullptr, box, FALSE);
    // check if atom numbers in Topology and Trajectory files match:
    natoms = read_first_x(oenv, &status, fnTRX, &t, &x, box);
    if (natoms > top.atoms.nr) { gmx_fatal(FARGS, "Topology (%d atoms) does not match trajectory (%d atoms)", top.atoms.nr, natoms); }
    // check if all atoms in selected index groups matching trajectory atoms:
    for(int ig=0; ig<ng; ig++) { if (isize[ig]>0) { check_index(nullptr, isize[ig], index[ig], nullptr, natoms); }}

    snew(qarr,isize[0]);
    for(int iq=0;iq<isize[0];iq++) qarr[iq]=Q4DEF;

    // initiate histogram arrays:
    binw = 1.0/(hist_bins-2);
    snew(sg_hist_tblk, hist_bins);
    for(i=0;i<hist_bins;i++)
    {
        sg_bin_centers[i] = (i-onehalf)*binw;
    }

    *nframes = 0;
    framenr  = 0;

    /* Loop over frames*/
    //TODO: add time average - accumulate x[] from NT read_next_x() and pass x_avg and box_avg to calc_tetra_order_hist()
    fprintf(stdout, "Processing frames:\n");
    do
    {
        // sample histogram every tblock steps (use the first frame without averaging for consistency)
        // Can calculate running average of the atom coordinates here before main calculator
        if (framenr % tblock == 0)
        {
            // Initialize histogram temporary storage:
            for(int iq=0;iq<isize[0];iq++) qarr[iq]=Q4DEF;
            (void)memset(sg_hist_tblk,0,hist_bins*sizeof(int)) ;

            // get q4 histogram on current frame:      
            calc_tetra_order_vals(top, pbcType, natoms, box, x, ng, isize, index, sflag, rlo, rhi, q4form, qarr);
            calc_hist( isize[0], qarr, hist_bins, binw, sg_hist_tblk);


            for (i = 0; i < hist_bins; i++)
            {
                sg_hist[i] += sg_hist_tblk[i] ;
            }

            // update collected frame count:
            (*nframes)++;

        }
        // update total frame count:
        framenr++;

    } while (read_next_x(oenv, status, &t, x, box));
    close_trx(status);

    sfree(sg_hist_tblk);
    sfree(qarr);
}
//-----------------------------------------------------------------------------------------------------------------------------
// Define groups and indices - forked from index.cpp: rd_groups, but uses predefined group indices  instead of querrying each group:
static void rd_groups_by_index(gmx::ArrayRef<const IndexGroup> indexGroups,
                               int                             ngrps,
                               int                             grpIdx[],
                               char*                           gnames[],
                               int                             isize[],
                               int*                            index[])
{

    for(int i=0; i<ngrps; i++) {
        if (grpIdx[i]>=0) {
            int gnr1 = grpIdx[i];
            gnames[i] = gmx_strdup(indexGroups[gnr1].name.c_str());
            isize[i]  = gmx::ssize(indexGroups[gnr1].particleIndices);
            snew(index[i], isize[i]);
            for (int j = 0; (j < isize[i]); j++)
            {
                index[i][j] = indexGroups[gnr1].particleIndices[j];
            }
        } else {
            isize[i] = 0; // predefined group was not found set size to 0 as a flag
        }
    }
}

// auxiliary to compute histogram - excluding the default -1.0 values:
static void calc_hist(  int     maxidx,
                        real*   qarr,
                        int     hist_bins,
                        real    hist_bin_width,
                        int*    sg_hist)
{
    int slindex_h;
    // Compute histogram
    for (int i = 0; (i < maxidx); i++) {
        if (qarr[i]>Q4DEF+0.1e-6) {
            if (qarr[i] < 0.0) {
                slindex_h = 0;
            } else if (qarr[i] > 1.0) {
                slindex_h = hist_bins-1;
            } else {
                slindex_h = static_cast<int>(std::floor(qarr[i] / hist_bin_width))+1;
            }
            sg_hist[slindex_h]++ ;
        }
    }
}

// Output of OTO values per frame:
static void write_q4vals(int frame, int natm_write, int*  index, real* qarr)
{
    FILE *ofile;
    char ofile_name[255];

    sprintf(ofile_name, "q4vals.%06d.dat",frame );
    fprintf(stdout, "Writing frame %06d data to %s\n",frame,ofile_name);

    ofile = gmx_ffopen(ofile_name, "w");
    for (int i = 0; i < natm_write; i++)
    {
        fprintf(ofile, "%d, %12.6f\n",index[i],qarr[i]);
    }
    gmx_ffclose(ofile);

}

// Output of hystogram:
static void writehist(int nframes,real* sg_bin_centers, int* sg_hist, int hist_bins, gmx::ArrayRef<const std::string> fnms)
{
    FILE *raw0;
    int   i;
    real  factor,bin_step;
    real  *sg_hist_prob, *sg_hist_tavg;

    if (nframes>0) {

        snew(sg_hist_tavg,hist_bins);
        snew(sg_hist_prob,hist_bins);

        // normalize by frame count:
        for (i = 0; i < hist_bins; i++)
        {
            sg_hist_tavg[i]=sg_hist[i]*1.0/nframes;
        }

        // normalization constant - only in the interval [0;1]:
        factor = 0.0;
        for (i = 1; i < hist_bins-1; i++)
        {
            factor+=sg_hist[i];
        }

        if (factor>0.000001) {
            for (i = 1; i < hist_bins-1; i++)
            {
                sg_hist_prob[i]=sg_hist[i]/factor;
            }
        }

        bin_step = sg_bin_centers[1]-sg_bin_centers[0];

        fprintf(stdout, "Writing histogram data to %s\n",fnms[0].c_str());
        raw0 = gmx_ffopen(fnms[0], "w");
        fprintf(raw0, "#bin_number\tbin_lower_edge\tbin_center\tbin_upper_edge\ttotal_counts\tcounts_averaged_by_frames\tbin_probability_[0:1]\n");
        for (i = 0; i < hist_bins; i++)
        {
            fprintf(raw0, "%10d\t%6.3f\t%6.3f\t%6.3f\t%12d\t%8.2f\t%12.6f\n",
                    i,
                    sg_bin_centers[i]-0.5*bin_step,
                    sg_bin_centers[i],
                    sg_bin_centers[i]+0.5*bin_step,
                    sg_hist[i],
                    sg_hist_tavg[i],
                    sg_hist_prob[i]);
        }
        gmx_ffclose(raw0);
        fprintf(stdout, "Finished histogrm output");

        sfree(sg_hist_tavg);
        sfree(sg_hist_prob);

    } else {
        fprintf(stdout,"Warning: No frames processed - empty output!\n");
    }
}
