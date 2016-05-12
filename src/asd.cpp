/*
asd - a program to quickly calculate pairwise individual allele sharing distances
    Copyright (C) 2011  Zachary A Szpiech (szpiechz@umich.edu)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <pthread.h>
//#include <cstdio>
#include <cstdlib>
#include <cctype>
#include "gzstream.h"
#include "param_t.h"
#include "asd-cli.h"
#include "asd-data.h"
#include "errlog.h"
#include "asd-dist.h"

using namespace std;


int main(int argc, char *argv[])
{
    param_t *params = getCLI(argc, argv);
    if (params == NULL) return 0;

    string outfile = params->getStringFlag(ARG_OUTFILE);
    LOG.init(outfile);
    LOG.log(getCommandLineString(argc, argv));
    LOG.log("Output file basename:", outfile);

    bool argerr = false;

    string filename = params->getStringFlag(ARG_FILENAME);
    string tped_filename = params->getStringFlag(ARG_TPED_FILENAME);
    string tfam_filename = params->getStringFlag(ARG_TFAM_FILENAME);
    bool STRU = (filename.compare(DEFAULT_FILENAME) != 0);
    bool TPED = (tped_filename.compare(DEFAULT_TPED_FILENAME) != 0);
    bool TFAM = (tfam_filename.compare(DEFAULT_TFAM_FILENAME) != 0);
    if (!check_file_type(STRU, TPED, TFAM)) {
        argerr = true;
        LOG.err("ERROR: Must specify either stru or tped/tfam.");
    }

    int nrows = params->getIntFlag(ARG_NROWS); //nchr
    if (!check_int_gt_0(nrows) && STRU) {
        argerr = true;
        LOG.err("ERROR: Number of chr must be > 0.");
    }
    LOG.log("nrows:", nrows);

    int ndrows = params->getIntFlag(ARG_NDROWS);
    if (!check_int_gt_0(ndrows) && STRU) {
        argerr = true;
        LOG.err("ERROR: Non-data rows must be > 0.");
    }
    LOG.log("ndrows:", ndrows);

    int ncols = params->getIntFlag(ARG_NCOLS);//nloci
    if (!check_int_gt_0(ncols) && STRU) {
        argerr = true;
        LOG.err("ERROR: Number of loci must be > 0.");
    }
    LOG.log("ncols:", ncols);

    int ndcols = params->getIntFlag(ARG_NDCOLS);
    if (!check_int_gt_0(ndcols) && STRU) {
        argerr = true;
        LOG.err("ERROR: Non-data columns must be > 0.");
    }
    LOG.log("ndcols:", ndcols);

    int sort = params->getIntFlag(ARG_SORT);
    if (!check_int_gt_0(sort) && STRU) {
        argerr = true;
        LOG.err("ERROR: Column to sort by must be > 0");
    }
    LOG.log("sort:", sort);

    if (!check_sort_ge_ndcols(sort, ndcols) && STRU) {
        argerr = true;
        LOG.err("ERROR: Must sort by a non-data column (sort >= ndcols).");
    }

    int num_threads = params->getIntFlag(ARG_THREAD);
    if (!check_int_gt_0(num_threads)) {
        argerr = true;
        LOG.err("ERROR: Must have a positive number of threads.");
    }
    LOG.log("threads:", num_threads);

    int nind = nrows / 2;
    bool PRINT_FULL = params->getBoolFlag(ARG_FULL);
    bool PRINT_FULL_LOG = params->getBoolFlag(ARG_FULL_LOG);
    if (!check_print_full(PRINT_FULL, PRINT_FULL_LOG)) {
        argerr = true;
        LOG.err("ERROR: Must choose only one of --full, --full-log.");
    }
    LOG.log("Print allele sharing distances:", PRINT_FULL || PRINT_FULL_LOG);
    LOG.log("Log transformed:", PRINT_FULL_LOG);

    bool CALC_ALL_IBS = params->getBoolFlag(ARG_CALC_IBS);
    LOG.log("Output IBS matricies:", CALC_ALL_IBS);

    bool CHECK_FILE = params->getBoolFlag(ARG_CHECK_FILE);
    LOG.log("Check STRU file:", CHECK_FILE);

    bool CHECK_FILE_DEEP = params->getBoolFlag(ARG_CHECK_FILE_DEEP);
    LOG.log("Check STRU file deep:", CHECK_FILE_DEEP);

    if (!check_file_check(CHECK_FILE, CHECK_FILE_DEEP, STRU)) {
        argerr = true;
        LOG.err("ERROR:", ARG_CHECK_FILE, false);
        LOG.err(" and", ARG_CHECK_FILE_DEEP, false);
        LOG.err(" are not supported with tped/tfam files.");
    }

    string STRU_MISSING = params->getStringFlag(ARG_STRU_MISSING);
    string TPED_MISSING = params->getStringFlag(ARG_TPED_MISSING);
    if (STRU) {
        LOG.log("STRU missing code:", STRU_MISSING);
    }
    else {
        LOG.log("TPED missing code:", TPED_MISSING);
    }
    //bool ASD = params->getBoolFlag(ARG_CALC_ASD);
    //bool FST = params->getBoolFlag(ARG_CALC_FST);

    if (argerr) return -1;

    bool FILE_STATUS_GOOD;
    if ((CHECK_FILE || CHECK_FILE_DEEP) && STRU)
    {
        FILE_STATUS_GOOD = checkFile(params);
        if (FILE_STATUS_GOOD) cerr << "File appears to be ok.\n";
        return -1;
    }

//    igzstream fin;
    structure_data *data;
    if (STRU)
    {
        //data = readData_ind_asd(fin, data, sort, ndcols, ndrows, nrows, ncols, STRU_MISSING);
        try {
            data = readData_stru(filename, sort, ndcols, ndrows, nrows, ncols, STRU_MISSING);
        }
        catch (...){
            return -1;
        }
    }
    else
    {
        nrows = 0;
        ncols = 0;
        /* if(ASD)*/
        try{
            //readData_ind_asd_tped_tfam(tped_filename, tfam_filename, data, nrows, ncols, TPED_MISSING);
        }
        catch (...)
        {
            return -1;
        }
        nind = nrows / 2;
    }

    init_storage(nind, CALC_ALL_IBS);

    unsigned int *NUM_PER_THREAD = make_thread_partition(num_threads, ncols);

    work_order_t *order;
    pthread_t *peer = new pthread_t[num_threads];
    unsigned int prev_index = 0;
    for (int i = 0; i < num_threads; i++)
    {
        order = new work_order_t;
        order->first_index = prev_index;
        order->last_index = prev_index + NUM_PER_THREAD[i];
        prev_index += NUM_PER_THREAD[i];
        order->stru_data = data;
        order->CALC_ALL_IBS = CALC_ALL_IBS;
        pthread_create(&(peer[i]),
                       NULL,
                       (void *(*)(void *))calc_pw_as_dist2,
                       (void *)order);

    }

    for (int i = 0; i < num_threads; i++) pthread_join(peer[i], NULL);

    finalize_calculations(nind, ncols, CALC_ALL_IBS);

    write_dist_matrix(outfile, nind, ncols, data->ind_names, PRINT_FULL, PRINT_FULL_LOG);

    if (CALC_ALL_IBS)
    {
        write_ibs_matrices(outfile, nind, ncols, data->ind_names, PRINT_FULL, PRINT_FULL_LOG);
    }

    delete [] NUM_PER_THREAD;
    delete [] peer;
    return 0;
}

