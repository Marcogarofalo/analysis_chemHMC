#define CONTROL
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// #include <time.h>
// #include <complex.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include "global.hpp"
#include "read.hpp"
#include "resampling_new.hpp"

#include "linear_fit.hpp"
#include "mutils.hpp"
#include "various_fits.hpp"

#include "correlators_analysis.hpp"
#include "functions_analysis_chemHMC.hpp"
struct kinematic kinematic_2pt;

generic_header read_head(FILE* stream) {
    int cs = 0;
    generic_header header;
    header.rs.resize(3);
    cs += fscanf(stream, "%lf\n", &header.rs[0]);//L
    cs += fscanf(stream, "%lf\n", &header.rs[1]);
    cs += fscanf(stream, "%lf\n", &header.rs[2]);
    cs += fscanf(stream, "%d\n", &header.size);// Nparticles
    header.mus.resize(1);
    cs += fscanf(stream, "%lf\n", &header.mus[2]);// mass
    cs += fscanf(stream, "%lf\n", &header.beta);// beta
    header.thetas.resize(5);
    cs += fscanf(stream, "%lf\n", &header.thetas[0]); //cutoff
    cs += fscanf(stream, "%lf\n", &header.thetas[1]); //eps
    cs += fscanf(stream, "%lf\n", &header.thetas[2]); //sigma
    cs += fscanf(stream, "%lf\n", &header.thetas[3]); //Lmax
    cs += fscanf(stream, "%lf\n", &header.thetas[4]); //bin_size

    header.T = (int)(header.thetas[3] / header.thetas[4]);
    header.ncorr = 1;
    cs += fscanf(stream, "%d\n", &header.Njack); //bin_size
    error(cs != 12, 1, "read_head", "error in reading the header");
    return header;
}
void write_header_g2(FILE* jack_file, generic_header head) {
    int fi = 0;
    fi += fwrite(&head.T, sizeof(int), 1, jack_file);
    fi += fwrite(&head.L, sizeof(int), 1, jack_file);
    int nmu = head.mus.size();
    fi += fwrite(&nmu, sizeof(int), 1, jack_file);
    for (double mu : head.mus) {
        fi += fwrite(&mu, sizeof(double), 1, jack_file);
    }
}

char** argv_to_options(char** argv) {
    char** option;
    option = (char**)malloc(sizeof(char*) * 7);
    option[0] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[1] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[2] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[3] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[4] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[5] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[6] = (char*)malloc(sizeof(char) * NAMESIZE);

    mysprintf(option[1], NAMESIZE, "read_plateaux"); // blind/see/read_plateaux
    mysprintf(option[2], NAMESIZE, "-p");            // -p
    mysprintf(option[3], NAMESIZE, argv[2]);         // path
    mysprintf(option[4], NAMESIZE, argv[6]);         // resampling
    mysprintf(option[5], NAMESIZE, "no");            // pdf
    mysprintf(option[6], NAMESIZE, argv[3]);         // infile
    return option;
}

void init_global_head(generic_header head) {
    file_head.l1 = head.L;
    file_head.l0 = head.T;
    file_head.l2 = head.L;
    file_head.l3 = head.L;
    file_head.nk = 2;
    file_head.musea = 0;
    file_head.k = (double*)malloc(sizeof(double) * file_head.nk * 2);
    file_head.k[0] = 0;
    file_head.k[1] = 0;
    file_head.k[2] = 0;
    file_head.k[3] = 0;

    file_head.nmoms = 1;
    file_head.mom = (double**)malloc(sizeof(double*) * file_head.nmoms);
    for (int i = 0; i < file_head.nmoms; i++) {
        file_head.mom[i] = (double*)malloc(sizeof(double) * 4);
        file_head.mom[i][0] = 0;
        file_head.mom[i][1] = 0;
        file_head.mom[i][2] = 0;
        file_head.mom[i][3] = 0;
    }
}

void read_twopt(FILE* stream, double*** to_write, generic_header head) {
    // write your function to read the data
    // int fi = 0;
    // for (int k = 0; k < head.ncorr; k++) {
    //     for (int t = 0; t < head.T;t++) {
    //         fi += fscanf(stream, "%lf  %lf", &to_write[k][t][0],
    //         &to_write[k][t][1]);
    //         // printf(" corr=%d t=%d %.12g   %.12g\n", k, t, to_write[k][t][0],
    //         to_write[k][t][1]);
    //     }
    // }
}
void how_many_confs_in_xyz(FILE* file, generic_header& head, int& confs) {
    int lines = 1;
    char c;
    int N;
    int tot = 0;
    tot += fscanf(file, "%d\n", &N);
    error(tot != 1, 1, "how_many_confs_in_xyz", "read: %d", tot);
    head.ncorr = N;
    /* count the newline characters */
    while ((c = fgetc(file)) != EOF) {
        if (c == '\n')
            lines++;
    }
    if (lines % (N + 2) != 0) {
        printf("error: input file contains %d lines\n", lines);
        printf("       the number of lines mus be a multiple of N+2=%d\n", N + 2);
        exit(1);
    }
    confs = lines / (N + 2);
    rewind(file);
}

void read_xyz(FILE* file, double**** to_write, generic_header head) {

    int count = 0;
    char id[1000];
    int N;
    char c;
    int tot = 0;
    for (int i = 0; i < head.Nconf; i++) {
        tot += fscanf(file, "%d\n", &N);
        error(N != i, 1, "read_RDF", "missmatch conf number read=%d expected=%d", N, i);
        for (int t = 0; t < head.T; t++) {
            tot += fscanf(file, "%lf   %lf \n",
                &to_write[i][0][t][1], &to_write[i][0][t][0]);
        }
    }
    error(tot != head.Nconf * (1 + 2 * head.T), 1, "how_many_confs_in_xyz",
        "read: %d", tot);
}
int main(int argc, char** argv) {
    error(argc != 7, 1, "nissa_mpcac ",
        "usage:./nissa_mpcac -p path file -bin $bin  jack/boot \n separate "
        "path and file please");

    char** option = argv_to_options(argv);

    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], option[6]);

    char namefile_plateaux[NAMESIZE];
    mysprintf(namefile_plateaux, NAMESIZE, "plateaux.txt");

    FILE* infile = open_file(namefile, "r");

    //////////////////////////////////// read and setup header
    generic_header head = read_head(infile);
    head.print_header();
    // head.print_header();
    init_global_head(head);

    //////////////////////////////////// setup jackboot and binning
    int confs = head.Njack;
    int bin = atoi(argv[5]);
    int Neff = confs / bin;
    int Njack;
    if (strcmp(argv[6], "jack") == 0) {
        Njack = Neff + 1;
        myres = new resampling_jack(Neff);
    }
    else if (strcmp(argv[6], "boot") == 0) {
        Njack = (Neff * 2 + 1);
        myres = new resampling_boot(Neff * 2);
    }
    else {
        Njack = 0;
        error(1 == 1, 1, "main", "argv[7]= %s is not jack or boot", argv[7]);
    }
    head.Nconf = head.Njack;
    head.Njack = Njack;
    //////////////////////////////////// setup output files
    mysprintf(namefile, NAMESIZE, "%s/out/%s_output", option[3], option[6]);
    printf("writing output in :\n %s \n", namefile);
    FILE* outfile = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/jackknife/%s_%s", option[3], option[4], option[6]);
    FILE* jack_file = open_file(namefile, "w+");
    head.write_header(jack_file);
    // write_header_g2(jack_file, head);

    //////////////////////////////////// confs
    double**** data = calloc_corr(confs, head.ncorr, head.T);

    printf("confs=%d\n", confs);
    printf("ncorr=%d\n", head.ncorr);
    // printf("kappa=%g\n", head.kappa);
    read_xyz(infile, data, head);
   
    double**** data_bin = binning(confs, head.ncorr, head.T, data, bin);
    double**** conf_jack = myres->create(Neff, head.ncorr, head.T, data_bin);
    free_corr(Neff, head.ncorr, head.T, data_bin);
    free_corr(confs, head.ncorr, head.T, data);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // print all the effective masses correlators
    // set the option to not read for a plateaux
    mysprintf(namefile, NAMESIZE, "%s/out/%s_meff_correlators", option[3],
        option[6]);
    FILE* outfile_meff_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_raw_correlators", option[3],
        option[6]);
    FILE* outfile_raw_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_shifted_correlators", option[3],
        option[6]);
    FILE* outfile_shifted_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_log_meff_shifted", option[3],
        option[6]);
    FILE* outfile_log_meff_shifted = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_gamma", option[3], option[6]);
    FILE* out_gamma = open_file(namefile, "w+");

    char save_option[NAMESIZE];
    sprintf(save_option, "%s", option[1]);
    sprintf(option[1], "blind");
    FILE* dev_null = open_file("/dev/null", "w");
    struct fit_type fit_info_silent;
    fit_info_silent.verbosity = -1;
    fit_info_silent.chi2_gap_jackboot = 1e+6;
    fit_info_silent.guess_per_jack = 0;

    for (int icorr = 0; icorr < head.ncorr; icorr++) {
        // log effective mass
        double* tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_meff_corr, icorr, "log", M_eff_log, dev_null,
            fit_info_silent);
        free(tmp_meff_corr);
        // raw correlator
        file_head.l0 = head.T*2;
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_raw_corr, icorr, "cor", identity, dev_null,
            fit_info_silent);
        free(tmp_meff_corr);
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_raw_corr, icorr, "cor", identity_im,
            dev_null, fit_info_silent);
        free(tmp_meff_corr);
        file_head.l0 = head.T;
        // shifted correlator
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_shifted_corr, icorr, "shift_cor", shift_corr,
            dev_null, fit_info_silent);
        free(tmp_meff_corr);
        // log_meff shifted correlator
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_log_meff_shifted, icorr, "log_shift",
            M_eff_log_shift, dev_null, fit_info_silent);
        free(tmp_meff_corr);
    }
    fit_info_silent.restore_default();
    sprintf(option[1], "%s", save_option); // restore option
    corr_counter = -1;


    ////////////////////////////////////////////////////////////
    // symmetrize
    ////////////////////////////////////////////////////////////
    // symmetrise_jackboot(Njack, 0, head.T, conf_jack);
    // symmetrise_jackboot(Njack, 1, head.T, conf_jack, -1);

    ////////////////////////////////////////////////////////////
    // start fitting
    //////////////////////////////

    // double* M_PS = plateau_correlator_function(
    //     option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
    //     namefile_plateaux, outfile, 0, "M_{PS}", identity, jack_file);
    // free(M_PS);
    // check_correlatro_counter(0);

    // eg of fit to correlator
    // struct fit_type fit_info;
    // struct fit_result fit_out;

    // fit_info.Nvar = 1;
    // fit_info.Npar = 1;
    // fit_info.N = 1;
    // fit_info.Njack = Njack;
    // fit_info.n_ext_P = 0;
    // // fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
    // // fit_info.ext_P[0] = something;
    // fit_info.function = constant_fit;
    // fit_info.linear_fit = true;
    // fit_info.T = head.T;
    // fit_info.corr_id = {};

    // // c++ 0 || r 1
    // struct fit_result fit_me_P5P5 = fit_fun_to_fun_of_corr(
    //     option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
    //     outfile, lhs_function_analysis_chemHMC_eg, "radial_distr", fit_info,
    //     jack_file);
    // check_correlatro_counter(0);
    // // free_fit_result(fit_info, fit_out);
    // fit_info.restore_default();

    ///////////// structure for fits
    //   data_all jackall;
    //   jackall.resampling = argv[6];
    //   jackall.ens = head.ncorr;
    //   jackall.en = new data_single[jackall.ens];
    //   for (int i = 0; i < jackall.ens; i++) {
    //     jackall.en[i].header = head;
    //     jackall.en[i].Nobs = 3;
    //     jackall.en[i].Njack = head.Njack;
    //     jackall.en[i].jack =
    //         (double **)malloc(sizeof(double *) * jackall.en[i].Nobs);
    //     for (int j = 0; j < Njack; j++) {
    //       jackall.en[i].jack[0][j] = conf_jack[j][i][0][0];
    //       jackall.en[i].jack[1][j] = conf_jack[j][i][1][0];
    //       jackall.en[i].jack[2][j] = conf_jack[j][i][2][0];
    //     }
    //   }
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // radial distribution function
    //   struct fit_type fit_info;
    //   fit_info.Npar = 3;
    //   fit_info.Nvar = 1;
    //   fit_info.Njack = head.Njack;
    //   // L=78.0838001  , max distance = 78.0838001*sqrt(3)=135.245109021
    //   fit_info.Prange={}
    //   fit_info.N = 136;
    //   fit_info.myen = std::vector<int>(jackall.ens);
    //   for (int n = 0; n < fit_info.myen.size(); n++)
    //     fit_info.myen[n] = n;
    //   fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() *
    //   fit_info.N,
    //                                fit_info.Njack);
    //   int count = 0;
    //   for (int n = 0; n < fit_info.N; n++) {
    //     for (int e : fit_info.myen) {
    //       for (int j = 0; j < fit_info.Njack; j++) {
    //         int ic = e % Nc;
    //         int is = e / Nc;
    //         fit_info.x[0][count][j] = 0; // mus1
    //         // fit_info.x[1][count][j] = head.mus[is + 1]; // mus1
    //       }
    //       count++;
    //     }
    //   }

    //   fit_info.corr_id = {0};
    //   fit_info.function = rhs_M_Ds_linear;
    //   fit_info.linear_fit = true;
    //   fit_info.covariancey = false;
    //   fit_info.verbosity = 0;
    //   mysprintf(namefit, NAMESIZE, "fit_MDs_vs_mu");

    //   fit_result fit_inter_MDs =
    //       fit_all_data(temp_argv, jackall_c, lhs_radial_distr_function,
    //       fit_info, namefit);
}