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
#include <cmath>
#include <cctype>
#include "param_t.h"
using namespace std;

const string ARG_CHECK_FILE = "--check-file";
const bool DEFAULT_CHECK_FILE = false;
const string HELP_CHECK_FILE = "Perform a quick check for file format integrity.  May miss some problems. (no calculations, no support for tped/tfam)";

const string ARG_CHECK_FILE_DEEP = "--check-deep";
const bool DEFAULT_CHECK_FILE_DEEP = false;
const string HELP_CHECK_FILE_DEEP = "Perform a longer check for file format integrity.  Includes all checks from --check-file. (no calculations, no support for tped/tfam)";

const string ARG_CALC_IBS = "--ibs-all";
const bool DEFAULT_CALC_IBS = false;
const string HELP_CALC_IBS = "Set to output all IBS sharing matrices.";

const string ARG_IBS_OUT = "--out";
const string DEFAULT_IBS_OUT = "outfile";
const string HELP_IBS_OUT = "Basename for output files.";

const string ARG_FILENAME = "--stru";
const string DEFAULT_FILENAME = "__none";
const string HELP_FILENAME = "The input data filename (stru format), data must be biallelic and coded 0/1 (-9 for missing data).";

const string ARG_TPED_FILENAME = "--tped";
const string DEFAULT_TPED_FILENAME = "__none";
const string HELP_TPED_FILENAME = "The input data filename (tped format), also requires a .tfam file (--tfam).  Data must be biallelic, but alleles can be arbitrarily coded.  Missing data must be -9.";

const string ARG_TFAM_FILENAME = "--tfam";
const string DEFAULT_TFAM_FILENAME = "__none";
const string HELP_TFAM_FILENAME = "The input data filename (tfam format), also requires a .tped file (--tped).";

const string ARG_NROWS = "--nchr";
const int DEFAULT_NROWS = 0;
const string HELP_NROWS = "Number of sampled chromosomes.";

const string ARG_NCOLS = "--nloci";
const int DEFAULT_NCOLS = 0;
const string HELP_NCOLS = "Number of loci.";

const string ARG_SORT = "--id-col";
const int DEFAULT_SORT = 1;
const string HELP_SORT = "Column where individual IDs are. (stru files only)";

const string ARG_NDCOLS = "--ndc";
const int DEFAULT_NDCOLS = 1;
const string HELP_NDCOLS = "Number of non-data columns. (stru files only)";

const string ARG_NDROWS = "--ndr";
const int DEFAULT_NDROWS = 1;
const string HELP_NDROWS = "Number of non-data rows. (stru files only)";

const string ARG_THREAD = "--threads";
const int DEFAULT_THREAD = 1;
const string HELP_THREAD = "Number of threads to spawn for faster calculation.";

const string ARG_FULL = "--full";
const bool DEFAULT_FULL = false;
const string HELP_FULL = "If set, prints the allele sharing matrix to stdout. Otherwise, outputs two r/2 x r/2 matrices. The first is the number of loci used for each pairwise comparison, and the second is the untransformed IBS matrix. In order to recover the allele sharing distance matrix divide the untransformed IBS counts by the number of loci in the appropriate cell, and then transform by 1-(ps)";

const string ARG_FULL_LOG = "--full-log";
const bool DEFAULT_FULL_LOG = false;
const string HELP_FULL_LOG = "Same as --full but transforms by -ln(ps).";

//const string ARG_MISSING = "-m";

const char DEL = ' ';
const string EMPTY_STRING = " ";

typedef struct
{
  short **data;
  string *locus_names;
  string *ind_names;
  //map<string,int> sample_size;
  int nloci;
  int nind;
} structure_data; 

typedef struct
{
  int first_index;
  int last_index;
  structure_data *stru_data;
  int missing;
  bool CALC_ALL_IBS;
} work_order_t;

typedef struct
{
  bool PRINT_FULL;
  bool PRINT_FULL_LOG;
  ostream *out;
  string *ind_names;
  int ncols;
  string type;
  int nind;
} output_order_t;

void output(void *order);
int search(string *s,int size,string key);
int put(string *s,int size,string key);

bool parse_cmd_line(int argc, char* argv[],
		    map<string,int> &argi,
		    map<string,string> &args,
		    map<string,bool> &argb);

void printHelp(map<string,int> &argi, map<string,string> &args,
	       map<string,bool> &argb);
void readData_ind_asd(ifstream &fin,structure_data &data,
	      int sort, int ndcols, int ndrows, int nrows, int ncols);
void readData_ind_asd_tped_tfam(ifstream &pedin, ifstream &famin, structure_data &data,
				int nind, int nloci);

short int* split_int(ifstream &fin, int fields);
string* split_str_str(int &size, const char *s, char c);
double proportion_shared(short A, short B);
void calc_pw_as_dist(void* work_order);

bool checkFile(param_t& params);
int countFields(string junk);
bool readData_Check(ifstream &fin,structure_data &data,
		    int sort, int ndcols, int ndrows, 
		    int nrows, int ncols);
void readData_pop_freq(ifstream &fin,structure_data &data,
		       int sort, int ndcols, int ndrows,
		       int nrows, int ncols);


double **DIST_MAT;
int **NUM_LOCI;
int **IBS_0_MAT;
int **IBS_1_MAT;
int **IBS_2_MAT;
pthread_mutex_t mutex_dist_mat = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_loci_mat = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_ibs_0 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_ibs_1 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_ibs_2 = PTHREAD_MUTEX_INITIALIZER;

int main(int argc, char* argv[])
{
  param_t params;
  params.addFlag(ARG_THREAD,DEFAULT_THREAD,"",HELP_THREAD);
  params.addFlag(ARG_CHECK_FILE,DEFAULT_CHECK_FILE,"",HELP_CHECK_FILE);
  params.addFlag(ARG_CHECK_FILE_DEEP,DEFAULT_CHECK_FILE_DEEP,"",HELP_CHECK_FILE_DEEP);
  params.addFlag(ARG_CALC_IBS,DEFAULT_CALC_IBS,"",HELP_CALC_IBS);
  params.addFlag(ARG_IBS_OUT,DEFAULT_IBS_OUT,"",HELP_IBS_OUT);
  params.addFlag(ARG_FILENAME,DEFAULT_FILENAME,"",HELP_FILENAME);
  params.addFlag(ARG_TPED_FILENAME,DEFAULT_TPED_FILENAME,"",HELP_TPED_FILENAME);
  params.addFlag(ARG_TFAM_FILENAME,DEFAULT_TFAM_FILENAME,"",HELP_TFAM_FILENAME);
  params.addFlag(ARG_NROWS,DEFAULT_NROWS,"",HELP_NROWS);
  params.addFlag(ARG_NCOLS,DEFAULT_NCOLS,"",HELP_NCOLS);
  params.addFlag(ARG_SORT,DEFAULT_SORT,"",HELP_SORT);
  params.addFlag(ARG_NDCOLS,DEFAULT_NDCOLS,"",HELP_NDCOLS);
  params.addFlag(ARG_NDROWS,DEFAULT_NDROWS,"",HELP_NDROWS);
  params.addFlag(ARG_FULL,DEFAULT_FULL,"",HELP_FULL);
  params.addFlag(ARG_FULL_LOG,DEFAULT_FULL_LOG,"",HELP_FULL_LOG);

  try
    {
      params.parseCommandLine(argc,argv);
    }
  catch (...)
    {
      return 1;
    }

  string outname = params.getStringFlag(ARG_IBS_OUT);
  string filename = params.getStringFlag(ARG_FILENAME);
  string tped_filename = params.getStringFlag(ARG_TPED_FILENAME);
  string tfam_filename = params.getStringFlag(ARG_TFAM_FILENAME);
  int nrows = params.getIntFlag(ARG_NROWS); //nchr
  int ndrows = params.getIntFlag(ARG_NDROWS);
  int ncols = params.getIntFlag(ARG_NCOLS);//nloci
  int ndcols = params.getIntFlag(ARG_NDCOLS);
  int sort = params.getIntFlag(ARG_SORT);
  int num_threads = params.getIntFlag(ARG_THREAD);
  int nind = nrows/2;
  bool PRINT_FULL = params.getBoolFlag(ARG_FULL);
  bool PRINT_FULL_LOG = params.getBoolFlag(ARG_FULL_LOG);
  bool quit = false;
  bool CALC_ALL_IBS = params.getBoolFlag(ARG_CALC_IBS);
  bool CHECK_FILE = params.getBoolFlag(ARG_CHECK_FILE);
  bool CHECK_FILE_DEEP = params.getBoolFlag(ARG_CHECK_FILE_DEEP);


  if(nrows <= 0)
    {
      cerr << "Number of rows must be > 0.\n";
      quit = true;
    }
  if(ncols <= 0)
    {
      cerr << "Number of cols must be > 0.\n";
      quit = true;
    }
  if(sort <= 0)
    {
      cerr << "Column to sort by must be > 0.\n";
      quit = true;
    }
  if(ndcols <= 0)
    {
      cerr << "Non-data columns must be > 0.\n";
      quit = true;
    }
  if(ndrows <= 0)
    {
      cerr << "Non-data rows must be > 0.\n";
      quit = true;
    }
  if(sort > ndcols)
    {
      cerr << "Must sort by a non-data column.\n";
      quit = true;
    }
  if(num_threads <= 0)
    {
      cerr << "Must have a positive number of threads.\n";
      quit = true;
    }
  if(num_threads > ncols)
    {
      cerr << "Number of threads must be < number of loci.\n";
      quit = true;
    }
  if(PRINT_FULL && PRINT_FULL_LOG)
    {
      cerr << "Must choose only one of --full, --full-log.\n";
      quit = true;
    }
  if(filename.compare(DEFAULT_FILENAME) == 0 && 
     tped_filename.compare(DEFAULT_TPED_FILENAME) == 0 && 
     tfam_filename.compare(DEFAULT_TFAM_FILENAME) == 0)
    {
      cerr << "Must specify a data file to read, either stru or tped/tfam.\n";
      quit = true;
    }

  if(filename.compare(DEFAULT_FILENAME) != 0 && 
     (tped_filename.compare(DEFAULT_TPED_FILENAME) != 0 || 
      tfam_filename.compare(DEFAULT_TFAM_FILENAME) != 0))
    {
      cerr << "Must specify only one type of data file to read, either stru or tped/tfam.\n";
      quit = true;
    }

  if(filename.compare(DEFAULT_FILENAME) == 0 && 
     ((tped_filename.compare(DEFAULT_TPED_FILENAME) == 0 && 
       tfam_filename.compare(DEFAULT_TFAM_FILENAME) != 0) ||
      (tped_filename.compare(DEFAULT_TPED_FILENAME) != 0 && 
       tfam_filename.compare(DEFAULT_TFAM_FILENAME) == 0)))
    {
      cerr << "Must specify both tped and tfam (or use stru).\n";
      quit = true;
    }
  if((CHECK_FILE || CHECK_FILE_DEEP) && (tped_filename.compare(DEFAULT_TPED_FILENAME) != 0 || 
					 tfam_filename.compare(DEFAULT_TFAM_FILENAME) != 0))
    {
      cerr << ARG_CHECK_FILE << " and " << ARG_CHECK_FILE_DEEP << " are not supported with tped/tfam files.\n";
      quit = true;
    }

  bool STRU_DATA = false;

  if(filename.compare(DEFAULT_FILENAME) != 0) STRU_DATA = true; 


  if(quit) return -1;

  bool FILE_STATUS_GOOD;
  if((CHECK_FILE || CHECK_FILE_DEEP) && (tped_filename.compare(DEFAULT_TPED_FILENAME) == 0 && 
      tfam_filename.compare(DEFAULT_TFAM_FILENAME) == 0))
    {
      FILE_STATUS_GOOD = checkFile(params);
      if(FILE_STATUS_GOOD) cerr << "File appears to be ok.\n";
      return -1;
    }  

  DIST_MAT = new double*[nind];
  for(int i = 0; i < nind;i++)
    {
      DIST_MAT[i] = new double[nind];
      for(int j = 0; j < nind;j++) DIST_MAT[i][j] = 0;
    }

  NUM_LOCI = new int*[nind];
  for(int i = 0; i < nind;i++)
    {
      NUM_LOCI[i] = new int[nind];
      for(int j = 0; j < nind;j++) NUM_LOCI[i][j] = 0;
    }

  if(CALC_ALL_IBS)
    {
      IBS_0_MAT = new int*[nind];
      for(int i = 0; i < nind;i++)
	{
	  IBS_0_MAT[i] = new int[nind];
	  for(int j = 0; j < nind;j++) IBS_0_MAT[i][j] = 0;
	}
      IBS_1_MAT = new int*[nind];
      for(int i = 0; i < nind;i++)
	{
	  IBS_1_MAT[i] = new int[nind];
	  for(int j = 0; j < nind;j++) IBS_1_MAT[i][j] = 0;
	}
      IBS_2_MAT = new int*[nind];
      for(int i = 0; i < nind;i++)
	{
	  IBS_2_MAT[i] = new int[nind];
	  for(int j = 0; j < nind;j++) IBS_2_MAT[i][j] = 0;
	}
    }


  ifstream fin,fin2;
  structure_data data;
  if(STRU_DATA)
    {
      fin.open(filename.c_str());
  
      if(fin.fail())
	{
	  cerr << "Could not open " << filename << " for reading.'n";
	  return -1;
	}
      readData_ind_asd(fin,data,sort,ndcols,ndrows,nrows,ncols);
    }
  else
    {
      fin.open(tped_filename.c_str());
      if(fin.fail())
	{
	  cerr << "Could not open " << tped_filename << " for reading.'n";
	  return -1;
	}

      fin2.open(tfam_filename.c_str());
      if(fin2.fail())
	{
	  cerr << "Could not open " << tfam_filename << " for reading.'n";
	  return -1;
	}      
      readData_ind_asd_tped_tfam(fin,fin2,data,nind,ncols);
    }
  

  data.nind = nrows/2;
  work_order_t *order;
  unsigned long int *NUM_PER_THREAD = new unsigned long int[num_threads];
  unsigned long int div = ncols/num_threads;

  for(int i = 0; i < num_threads; i++)
    {
      NUM_PER_THREAD[i] = 0;
      NUM_PER_THREAD[i] += div;
    }

  for(int i = 0; i < ncols%num_threads;i++)
    {
      NUM_PER_THREAD[i]++;
    }
 
  pthread_t *peer = new pthread_t[num_threads];
  unsigned long int prev_index = 0;
  for(int i = 0; i < num_threads; i++)
    {
      order = new work_order_t;
      order->first_index = prev_index;
      order->last_index = prev_index+NUM_PER_THREAD[i];
      prev_index += NUM_PER_THREAD[i];
      order->stru_data = &data;
      order->CALC_ALL_IBS = CALC_ALL_IBS;
      pthread_create(&(peer[i]),
		     NULL,
		     (void *(*)(void*))calc_pw_as_dist,
		     (void *)order);
      
    }


  for(int i = 0; i < num_threads; i++)
    {
      pthread_join(peer[i],NULL);
    }


  for(int i = 0; i < nind;i++)
    {
      for(int j = i; j < nind;j++)
	{
	  if (i==j)
	    {
	      DIST_MAT[i][j] = ncols+NUM_LOCI[i][j];
	      //NUM_LOCI[i][j] = ncols+NUM_LOCI[i][j];
	      if(CALC_ALL_IBS)
		{
		  IBS_0_MAT[i][j] = 0;
		  IBS_1_MAT[i][j] = 0;
		  IBS_2_MAT[i][j] = ncols+NUM_LOCI[i][j];
		}
	    }
	  else
	    {
	      DIST_MAT[j][i] = DIST_MAT[i][j];
	      NUM_LOCI[j][i] = NUM_LOCI[i][j];
	      if(CALC_ALL_IBS)
		{
		  IBS_0_MAT[j][i] = IBS_0_MAT[i][j];
		  IBS_1_MAT[j][i] = IBS_1_MAT[i][j];
		  IBS_2_MAT[j][i] = IBS_2_MAT[i][j];
		}
	    }
	}
    }


  //I should have checked if streams were threadsafe before I bothered with this
  //But they aren't, or I missed something, so I join threads immediately after
  //creating

  ostream **out = new ostream*[4];
  ofstream outfile_dist, outfile_ibs0, outfile_ibs1, outfile_ibs2;
  if(outname.compare("outfile") == 0)
    {
      out[0] = &cout;
      if(CALC_ALL_IBS)
	{
	  out[1] = &cout;
	  out[2] = &cout;
	  out[3] = &cout;
	}
    }
  else
    {
      string dist_fname = outname+".dist";
      outfile_dist.open(dist_fname.c_str());
      out[0] = &outfile_dist;
      if(CALC_ALL_IBS)
	{
	  string ibs0_fname = outname+".ibs0";
	  string ibs1_fname = outname+".ibs1";
	  string ibs2_fname = outname+".ibs2";
	  outfile_ibs0.open(ibs0_fname.c_str());
	  out[1] = &outfile_ibs0;
	  outfile_ibs1.open(ibs1_fname.c_str());
	  out[2] = &outfile_ibs1;
	  outfile_ibs2.open(ibs2_fname.c_str());
	  out[3] = &outfile_ibs2;
	}
    }

  pthread_t output_peer;
  output_order_t **output_order = new output_order_t*[4];
  string types[4] = {"dist","ibs0","ibs1","ibs2"};
 

  for(int i = 0; i < 4; i ++)
    {
      output_order[i] = new output_order_t;
      output_order[i]->PRINT_FULL = PRINT_FULL;
      output_order[i]->PRINT_FULL_LOG = PRINT_FULL_LOG;
      output_order[i]->ind_names = data.ind_names;
      output_order[i]->ncols = ncols;
      output_order[i]->nind = nind;
    }
 
  
  for(int i = 0; i < 4; i++)
    {
      if(i > 0 && !CALC_ALL_IBS) continue;
      output_order[i]->out = out[i];
      output_order[i]->type = types[i];
      pthread_create(&(output_peer),
		     NULL,
		     (void *(*)(void*))output,
		     (void *)output_order[i]);
      //join after create because this was all useless
      pthread_join(output_peer,NULL);
    }
 

  delete [] NUM_PER_THREAD;
  delete [] peer;
  delete [] data.locus_names;

  for(int i = 0; i < 4; i++)
    {
      delete output_order[i];
    }
  delete [] output_order;

  return 0;
}


bool readData_Check(ifstream &fin,structure_data &data,
	      int sort, int ndcols, int ndrows, 
	      int nrows, int ncols)
{
  bool FILE_STATUS = true;

  string line;
  int nind = nrows/2;
  data.nind = nind;
  data.data = new short*[nind];
  data.ind_names = new string[nind];

  for (int i = 0; i < nind; i++)
    {
      data.ind_names[i] = EMPTY_STRING;
    }

  getline(fin,line);
  int size;
  data.locus_names = split_str_str(size,line.c_str(),DEL);

  size = ncols;
  data.nloci = size;

  for(int i = 1; i < ndrows;i++)
    {
      getline(fin,line);
    }

  string key;
  string field;
  short int *tmp;
  //int **block;
  //double *block;
  short tmp_dbl;
  int index;


  short *ind_count = new short[nind];

  for(int i = 0; i < nind; i++) ind_count[i] = 0;

  for(int row = 0; row < nrows;row++)
    {
      tmp = NULL;
      //block = NULL;
      for(int i = 1; i <= ndcols; i++)
	{
	  fin >> field;
	  if(i == sort) key = field;
	}

      index = search(data.ind_names,nind,key);
      tmp = split_int(fin,ncols);

      if(index >= 0)
	{
	  ind_count[index]++;

	  if(ind_count[index] > 2)
	    {
	       FILE_STATUS &= false;
	       cerr << "Individual " << key << " seen "
		    << ind_count[index] << "times.\n";
	    }

	  for(int i = 0; i < size; i++)
	    {

	      if(tmp[i] != 1 && tmp[i] != 0 && tmp[i] != -9)
		{
		  FILE_STATUS &= false;
		  cerr << "LINE " << row + ndrows << " COL "
		       << i+ndcols << ": Allele '" << tmp[i] 
		       << "' is not 0/1/-9.\n";
		}
	    }
	}
      else
	{
	  index = put(data.ind_names,nind,key);
	  if(index > nind-1)
	    {
	      cerr << "Found more than " << nind << " individuals.";
	      exit(-1);
	    }
	  ind_count[index]++;

	  for(int i = 0; i < size; i++)
	    {
	      if(tmp[i] != 1 && tmp[i] != 0 && tmp[i] != -9)
		{
		  FILE_STATUS &= false;
		  cerr << "LINE " << row + ndrows << " COL "
		       << i << ": Allele " << tmp[i] << " is not 0/1/-9.\n";
		}
	    }
	}
      delete [] tmp;
    }

  return FILE_STATUS;
}

bool checkFile(param_t& params)
{

  string outname = params.getStringFlag(ARG_IBS_OUT);
  string filename = params.getStringFlag(ARG_FILENAME);
  int nrows = params.getIntFlag(ARG_NROWS);
  int ndrows = params.getIntFlag(ARG_NDROWS);
  int ncols = params.getIntFlag(ARG_NCOLS);
  int ndcols = params.getIntFlag(ARG_NDCOLS);
  int sort = params.getIntFlag(ARG_SORT);
  int num_threads = params.getIntFlag(ARG_THREAD);
  int nind = nrows/2;
  bool PRINT_FULL = params.getBoolFlag(ARG_FULL);
  bool PRINT_FULL_LOG = params.getBoolFlag(ARG_FULL_LOG);
  bool quit = false;
  bool CALC_ALL_IBS = params.getBoolFlag(ARG_CALC_IBS);
  bool CHECK_FILE = params.getBoolFlag(ARG_CHECK_FILE);
  bool CHECK_FILE_DEEP = params.getBoolFlag(ARG_CHECK_FILE_DEEP);
  bool FILE_STATUS = true;

  ifstream fin;
  fin.open(filename.c_str());
  
  if(fin.fail())
    {
      cerr << "Could not open " << filename << " for reading.'n";
      return 0;
    }
  
  string junk;
  int fields;
  int counter = 1;
  int ndr_counter = ndrows;
  int obs_cols;
  int obs_rows = 0;
  do
    {
      getline(fin,junk);
      fields = countFields(junk);
      if(fields == 0 && counter < nrows+ndrows)
	{
	  FILE_STATUS &= false;
	  cerr << "LINE " << counter << ": "
	       << "Blank line found.\n";
	}
      else if(fields == 0)
	{
	}
      else if(ndr_counter > 0)
	{
	  if(counter == 1)
	    {
	      obs_cols = fields;
	      if(obs_cols != ncols)
		{
		  FILE_STATUS &= false;
		  cerr << "LINE " << counter << ": ";
		  cerr << "Found " << obs_cols << " loci. Expected " 
		       << ncols << ".\n";;
		}
	    }
	  else
	    {
	      if(fields != ncols)
		{
		  FILE_STATUS &= false;
		  cerr << "LINE " << counter << ": ";
		  cerr << "Found " << fields << " fields in an assumed "
		       << "header row. Expected " << ncols << ".\n";;
		}
	    }
	  ndr_counter--;
	}
      else
	{
	  if(fields-ncols != ndcols)
	    {
	      FILE_STATUS &= false;
	      cerr << "LINE " << counter << ": ";
	      cerr << "Found " << fields << " fields. Expected " 
		   << ndcols << " headers + " 
		   << ncols << " data = " 
		   << ndcols+ncols << " fields.\n";
	    }
	  obs_rows++;
	}
      
      counter++;
    } while(fin.good());

  if(obs_rows != nrows)
    {
      FILE_STATUS &= false;
      cerr << "Found " << obs_rows << " lines of data. "
	   << "Expected " << nrows << ".\n";
    }

  fin.close();

  //Instead of reopening file, we can just seekg() back to the beginning.
  //Change this at some point.

  short** data;
  
  if(FILE_STATUS && CHECK_FILE_DEEP)
    {
      fin.open(filename.c_str());
      structure_data data;
      FILE_STATUS &= readData_Check(fin,data,sort,ndcols,ndrows,nrows,ncols);
      fin.close();
    }


  return FILE_STATUS;
}


int countFields(string str)
{
  string::iterator it;
  bool inword = false;
  bool space = false;
  int counter = 0;
  for ( it=str.begin() ; it < str.end(); it++ )
    {
      space = isspace(*it);
      if(!inword && !space)
	{
	  inword = true;
	  counter++;
	}
      else if(inword && space)
	{
	  inword = false;
	}
    }

  return counter;
}


//Yes, this is really dumb for a lot of reasons.
void output(void *order)
{
  output_order_t *p = (output_order_t*)order;
  ostream *out = p->out;
  bool PRINT_FULL = p->PRINT_FULL;
  bool PRINT_FULL_LOG = p->PRINT_FULL_LOG;
  string *ind_names = p->ind_names;
  int ncols = p->ncols;
  string type = p->type;
  int nind = p->nind;

  if(type.compare("dist") == 0)
    {
      if(!PRINT_FULL && !PRINT_FULL_LOG)
	{
	  for(int i = 0; i < nind;i++)
	    {
	      for(int j = 0; j < nind;j++)
		{
		  *out << double(ncols)+double(NUM_LOCI[i][j]) << " ";   
		}
	      *out << endl;
	    }
	  *out << endl;
	}
      
      for (int i = 0; i < nind; i++)
	{
	  *out << ind_names[i] << " ";
	}
      *out << endl;
      
      for(int i = 0; i < nind; i++)
	{
	  *out << ind_names[i] << " ";
	  for (int j = 0; j < nind ; j++)
	    {
	      if(PRINT_FULL_LOG)
		{
		  *out << 0-log(double(DIST_MAT[i][j]) / 
				     (double(ncols)+double(NUM_LOCI[i][j]))) 
			    << " ";
		}
	      else if(PRINT_FULL)
		{
		  *out << 1-(double(DIST_MAT[i][j]) / 
				  (double(ncols)+double(NUM_LOCI[i][j]))) 
			    << " ";
		}		  
	      else *out << DIST_MAT[i][j] << " ";
	    }
	  *out << endl;
	}
    }
  else if(type.compare("ibs0") == 0)
    {
      if(!PRINT_FULL && !PRINT_FULL_LOG)
	{
	  for(int i = 0; i < nind;i++)
	    {
	      for(int j = 0; j < nind;j++)
		{
		  *out << double(ncols)+double(NUM_LOCI[i][j]) << " ";   
		}
	      *out << endl;
	    }
	  *out << endl;
	}
      
      for (int i = 0; i < nind; i++)
	{
	  *out << ind_names[i] << " ";
	}
      *out << endl;
      
      for(int i = 0; i < nind; i++)
	{
	  *out << ind_names[i] << " ";
	  for (int j = 0; j < nind ; j++)
	    {
	      if(PRINT_FULL_LOG || PRINT_FULL)
		{
		  *out << double(IBS_0_MAT[i][j])/
		    (double(ncols)+double(NUM_LOCI[i][j]))<< " ";
		}
	      else *out << IBS_0_MAT[i][j] << " ";

	    }
	  *out << endl;
	}
    }
  else if(type.compare("ibs1") == 0)
    {
      if(!PRINT_FULL && !PRINT_FULL_LOG)
	{
	  for(int i = 0; i < nind;i++)
	    {
	      for(int j = 0; j < nind;j++)
		{
		  *out << double(ncols)+double(NUM_LOCI[i][j]) << " ";   
		}
	      *out << endl;
	    }
	  *out << endl;
	}
      
      for (int i = 0; i < nind; i++)
	{
	  *out << ind_names[i] << " ";
	}
      *out << endl;
      
      for(int i = 0; i < nind; i++)
	{
	  *out << ind_names[i] << " ";
	  for (int j = 0; j < nind ; j++)
	    {
	      if(PRINT_FULL_LOG || PRINT_FULL)
		{
		  *out << double(IBS_1_MAT[i][j])/
		    (double(ncols)+double(NUM_LOCI[i][j]))<< " ";
		}
	      else *out << IBS_1_MAT[i][j] << " ";
	      
	    }
	  *out << endl;
	}
    }
  else if(type.compare("ibs2") == 0)
    {
      if(!PRINT_FULL && !PRINT_FULL_LOG)
	{
	  for(int i = 0; i < nind;i++)
	    {
	      for(int j = 0; j < nind;j++)
		{
		  *out << double(ncols)+double(NUM_LOCI[i][j]) << " ";   
		}
	      *out << endl;
	    }
	  *out << endl;
	}
      
      for (int i = 0; i < nind; i++)
	{
	  *out << ind_names[i] << " ";
	}
      *out << endl;
      
      for(int i = 0; i < nind; i++)
	{
	  *out << ind_names[i] << " ";
	  for (int j = 0; j < nind ; j++)
	    {
	      if(PRINT_FULL_LOG || PRINT_FULL)
		{
		  *out << double(IBS_2_MAT[i][j])/
		    (double(ncols)+double(NUM_LOCI[i][j]))<< " ";
		}
	      else *out << IBS_2_MAT[i][j] << " ";
	      
	    }
	  *out << endl;
	}
    }
  else
    {
      cerr << "UNKNOWN TYPE: " << type << endl;
      exit(-1);
    }
  return;
}

void calc_pw_as_dist(void* order)
{
  work_order_t *p = (work_order_t*)order;
  //map<string,double*>::iterator key;
  //map<string,double*> *data = p->stru_data->data;
  short **data = p->stru_data->data;
  int size = p->stru_data->nind;
  //string *key_list = new string[size];
  //int i = 0;
  /*
  for (key = data->begin(); key != data->end(); key++)
    {
      key_list[i] = (*key).first;
      i++;
    }
  */
  short A, B;

  double ps;
  double* row = NULL;
  int* num_loci = NULL;
  int* ibs0 = NULL;
  int* ibs1 = NULL;
  int* ibs2 = NULL;
  for(int j = 0; j < size; j++)
    {
      row = new double[size-j];
      num_loci = new int[size-j];
      if(p->CALC_ALL_IBS)
	{
	  ibs0 = new int[size-j];
	  ibs1 = new int[size-j];
	  ibs2 = new int[size-j];
	}
      for(int k = j; k < size; k++)
	{
	  row[k-j]=0;
	  num_loci[k-j]=0;
	  if(p->CALC_ALL_IBS)
	    {
	      ibs0[k-j]=0;
	      ibs1[k-j]=0;
	      ibs2[k-j]=0;
	    }
	  for(int l = p->first_index; l < p->last_index;l++)
	    {
	      if(j==k)
		{
		  A = data[j][l];
		  //B = data[k][l];
		  if(A < 0 /*|| B < 0*/)
		    {
		      num_loci[k-j]--;
		    }
		}
	      else
		{
		  A = data[j][l];
		  B = data[k][l];
		  if(A < 0 || B < 0)
		    {
		      num_loci[k-j]--;
		    }
		  else
		    {
		      ps = proportion_shared(A,B);
		      row[k-j]+=ps;
		      if(p->CALC_ALL_IBS)
			{
			  if(ps == 1) ibs2[k-j]++;
			  if(ps == 0.5) ibs1[k-j]++;
			  if(ps == 0) ibs0[k-j]++;
			}
		    }
		}
	    }
	}
     
      pthread_mutex_lock(&mutex_dist_mat);
      for(int m = j; m < size; m++)  DIST_MAT[j][m]+=double(row[m-j]);
      pthread_mutex_unlock(&mutex_dist_mat);
      
      pthread_mutex_lock(&mutex_loci_mat);
      for(int m = j; m < size; m++)  NUM_LOCI[j][m]+=num_loci[m-j];
      pthread_mutex_unlock(&mutex_loci_mat);

      if(p->CALC_ALL_IBS)
	{
	  pthread_mutex_lock(&mutex_ibs_0);
	  for(int m = j; m < size; m++)  IBS_0_MAT[j][m]+=ibs0[m-j];
	  pthread_mutex_unlock(&mutex_ibs_0);
	  
	  pthread_mutex_lock(&mutex_ibs_1);
	  for(int m = j; m < size; m++)  IBS_1_MAT[j][m]+=ibs1[m-j];
	  pthread_mutex_unlock(&mutex_ibs_1);
	  
	  pthread_mutex_lock(&mutex_ibs_2);
	  for(int m = j; m < size; m++)  IBS_2_MAT[j][m]+=ibs2[m-j];
	  pthread_mutex_unlock(&mutex_ibs_2);
	}

      delete [] num_loci;
      delete [] row;
      if(p->CALC_ALL_IBS)
	{
	  delete [] ibs0;
	  delete [] ibs1;
	  delete [] ibs2;
	}
    }

  //delete [] key_list;
  delete p;
  return;

}

double proportion_shared(short A, short B)
{
  if(abs(A - B) == 0) return 1;
  if(abs(A - B) == 1) return 0.5;
  return 0;
}

void readData_ind_asd(ifstream &fin,structure_data &data,
	      int sort, int ndcols, int ndrows, 
	      int nrows, int ncols)
{
  string line;
  int nind = nrows/2;
  data.nind = nind;
  data.data = new short*[nind];
  data.ind_names = new string[nind];

  for (int i = 0; i < nind; i++)
    {
      data.ind_names[i] = EMPTY_STRING;
    }

  getline(fin,line);
  int size;
  data.locus_names = split_str_str(size,line.c_str(),DEL);

  size = ncols;
  data.nloci = size;

  for(int i = 1; i < ndrows;i++)
    {
      getline(fin,line);
    }

  string key;
  string field;
  short int *tmp;
  //int **block;
  //double *block;
  short tmp_dbl;
  int index;
  for(int row = 0; row < nrows;row++)
    {
      tmp = NULL;
      //block = NULL;
      for(int i = 1; i <= ndcols; i++)
	{
	  fin >> field;
	  if(i == sort) key = field;
	}

      index = search(data.ind_names,nind,key);

      //getline(fin,line);
      //cerr << line << endl;
      tmp = split_int(fin,ncols);

      if(index >= 0)
	{
	  for(int i = 0; i < size; i++)
	    {
	      tmp_dbl = data.data[index][i];
	      //(*data.data)[key][1][i] = tmp[i];
	      if(tmp_dbl >= 0 && tmp[i] >= 0)
		{
		  data.data[index][i] += tmp[i];
		}
	      else
		{
		  data.data[index][i] = -9;
		}
	    }
	}
      else
	{
	  index = put(data.ind_names,nind,key);
	  data.data[index] = new short[size];
	  //block = new double[size];
	  /*
	  block = new int[2];
	  block[0] = new int[size];
	  block[1] = new int[size];
	  */	  
	  for(int i = 0; i < size; i++)
	    {
	      data.data[index][i] = tmp[i];
	    }
	  //data.data[index] = block;
	}
      delete [] tmp;
    }

  return;
}

void readData_ind_asd_tped_tfam(ifstream &pedin, ifstream &famin, structure_data &data,
				int nind, int nloci)
{
  string line;
  string junk;
  //int nind = nrows/2;
  
  data.nind = nind;
  data.data = new short*[nind];
  for(int i = 0; i < nind; i++) data.data[i] = new short[nloci];
  data.ind_names = new string[nind];
  
  for(int i = 0; i < nind; i++)
    {
      famin >> junk;
      famin >> data.ind_names[i];
      getline(famin,line);
    }
  
  data.nloci = nloci;
  data.locus_names = new string[nloci];
  
  string *zeroAllele = new string[nloci];
  for (int i = 0; i < nloci; i++) zeroAllele[i] = "-9";
  string allele1, allele2;
  short alleleCount = 0;
  for (int locus = 0; locus < nloci; locus++)
    {
      pedin >> junk;
      pedin >> data.locus_names[locus];
      pedin >> junk;
      pedin >> junk;

      for(int ind = 0; ind < nind; ind++)
	{
	  pedin >> allele1;
	  pedin >> allele2;
	  if(allele1.compare("-9") == 0 || allele2.compare("-9") == 0)
	    {
	      data.data[ind][locus] = -9;
	    }
	  else if(zeroAllele[locus].compare("-9") == 0)
	    {
	      alleleCount = 0;
	      zeroAllele[locus] = allele1;
	      if(allele2.compare(zeroAllele[locus]) != 0) alleleCount++;
	      data.data[ind][locus] = alleleCount;
	    }
	  else
	    {
	      alleleCount = 0;
	      if(allele1.compare(zeroAllele[locus]) != 0) alleleCount++;
	      if(allele2.compare(zeroAllele[locus]) != 0) alleleCount++;
	      data.data[ind][locus] = alleleCount;
	    }
	}
    }

  delete [] zeroAllele;

  return;
}


void readData_pop_freq(ifstream &fin,structure_data &data,
	      int sort, int ndcols, int ndrows, 
	      int nrows, int ncols)
{
  string line;
  int nind = nrows/2;
  data.nind = nind;
  data.data = new short*[nind];
  data.ind_names = new string[nind];

  for (int i = 0; i < nind; i++)
    {
      data.ind_names[i] = EMPTY_STRING;
    }

  getline(fin,line);
  int size;
  data.locus_names = split_str_str(size,line.c_str(),DEL);

  size = ncols;
  data.nloci = size;

  for(int i = 1; i < ndrows;i++)
    {
      getline(fin,line);
    }

  string key;
  string field;
  short int *tmp;
  //int **block;
  //double *block;
  short tmp_dbl;
  int index;
  for(int row = 0; row < nrows;row++)
    {
      tmp = NULL;
      //block = NULL;
      for(int i = 1; i <= ndcols; i++)
	{
	  fin >> field;
	  if(i == sort) key = field;
	}

      index = search(data.ind_names,nind,key);

      //getline(fin,line);
      //cerr << line << endl;
      tmp = split_int(fin,ncols);

      if(index >= 0)
	{
	  for(int i = 0; i < size; i++)
	    {
	      tmp_dbl = data.data[index][i];
	      //(*data.data)[key][1][i] = tmp[i];
	      if(tmp_dbl >= 0 && tmp[i] >= 0)
		{
		  data.data[index][i] += tmp[i];
		}
	      else
		{
		  data.data[index][i] = -9;
		}
	    }
	}
      else
	{
	  index = put(data.ind_names,nind,key);
	  data.data[index] = new short[size];
	  //block = new double[size];
	  /*
	  block = new int[2];
	  block[0] = new int[size];
	  block[1] = new int[size];
	  */	  
	  for(int i = 0; i < size; i++)
	    {
	      data.data[index][i] = tmp[i];
	    }
	  //data.data[index] = block;
	}
      delete [] tmp;
    }

  return;
}


int search(string *s,int size,string key)
{
  for (int i = 0; i < size; i++)
    {
      if(s[i].compare(key) == 0) return i;
    }
  return -9;
}

int put(string *s,int size,string key)
{
  for(int i = 0; i < size; i ++)
    {
      if(s[i].compare(EMPTY_STRING) == 0)
	{
	  s[i] = key;
	  return i;
	}
    }
  cerr << "ERROR: This shouldn't be possible.  Try --check-file or --check-deep.\n";
  exit(0);
}

// Split a NULL-terminated character array on a given character into
// a vector of strings
// The vector is passed by reference and cleared each time
// The number of strings split out is returned
short int* split_int(ifstream &fin, int fields)
{
  /*
  vector<string> v;
  while (true)
    {
      const char* begin = s;
      
      while (*s != c && *s) { ++s; }
      
      v.push_back(string(begin, s));
      
      if (!*s)
	{
	  break;
	}
      
      if (!*++s)
	{
	  //v.push_back("");
	  break;
	}
    }
  */
  short int *a = new short int[fields];
  for (int i = 0; i < fields; i++)
    {
      fin >> a[i];
      //cout << a[i] << endl;
    }
  /*
  for(int i = 0; i < v.size();i++)
    {
      a[i]=atoi(v[i].c_str());
    }
  */
  return a;
}

string* split_str_str(int &size, const char* s, char c)
{
  vector<string> v;
  while (true)
    {
      const char* begin = s;
      
      while (*s != c && *s) { ++s; }
      
      v.push_back(string(begin, s));
      
      if (!*s)
	{
	  break;
	}
      
      if (!*++s)
	{
	  v.push_back("");
	  break;
	}
    }
  string * x = new string[v.size()];
  
  for(int i = 0; i < v.size();i++)
    {
      x[i]=v[i];
    } 

  size = v.size();
  return x;
}

/*
void printHelp(map<string,int> &argi, map<string,string> &args,
	       map<string,bool> &argb)
{
  map<string,int>::iterator i1;
  map<string,string>::iterator i2;
  map<string,bool>::iterator i3;
  cout << "\n\n++++++++++Command Line Arguments++++++++++\n\n";

  for (i1=argi.begin();i1!=argi.end();i1++)
    {
      cout << (*i1).first<< " int\n";
    }
   for (i2=args.begin();i2!=args.end();i2++)
    {
      cout << (*i2).first<< " string\n";
    }
   for (i3=argb.begin();i3!=argb.end();i3++)
    {
      cout << (*i3).first<< " bool\n";
    }
  
   return;
}


bool parse_cmd_line(int argc, char* argv[],
		    map<string,int> &argi, 
		    map<string,string> &args,
		    map<string,bool> &argb)
{

  argb[ARG_CHECK_FILE_DEEP] = false;
  argb[ARG_CHECK_FILE] = false;
  args[ARG_IBS_OUT] = "outfile";
  args[ARG_FILENAME] = "infile";
  argi[ARG_NROWS] = 0;
  argi[ARG_NCOLS] = 0;
  argi[ARG_SORT] = 1;
  argi[ARG_NDCOLS] = 0;
  argi[ARG_NDROWS] = 0;
  argi[ARG_THREAD] = 1;
  argb[ARG_FULL] = false;
  argb[ARG_FULL_LOG] = false;
  argb[ARG_CALC_IBS] = false;
  
  
  for (int i = 1; i < argc;i++)
    {
      if(argi.count(argv[i]) > 0)
	{
	  argi[argv[i]] = atoi(argv[i+1]);
	}
      else if(args.count(argv[i]) > 0)
	{
	  args[argv[i]] = argv[i+1];
	}
      else if(argb.count(argv[i]) > 0)
	{
	  argb[argv[i]] = true;
	}
      else if(ARG_HELP.compare(argv[i]) == 0)
	{
	  return 1;
	}
    }

  return 0;
}

*/
