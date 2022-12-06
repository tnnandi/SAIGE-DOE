#include <iostream>
#include <fstream>
#include <vector>
#include <mpi.h>
#include <math.h>


//proof of concept code for parallel file writing on Summit (tested on Polaris)
// For some reason, the code throws error (after writing out the files) when writing in parallel using 4 process (but has no issues for other number of processes)

/*
qsub -I -A covid-ct -l select=1:system=polaris -l filesystems=home:eagle -l walltime=60:00 -q debug
module load mpiwrappers/cray-mpich-llvm
mpic++ <file_name>.cpp
mpirun -np <n_proc> ./a.out
*/


int main(int argc, char **argv){

        int *buffer; // buffer (local to a process) to hold the file content
        int rank;
        int size;
        int error;
        int mpi_chunk_size; // number of rows in each MPI process
        int mpi_k_start;    // starting row index in each process  
        int mpi_k_end;      // end row index in each process

	// some memory issue exits the porgram for the write_file_in_parallel case but does the job. So, for the other options turn write_file_in_parallel off
	int write_files_from_each_process = 0;
	int write_file_from_root          = 0;
	int write_file_in_parallel        = 1;

        MPI_File fh;
        MPI_Status status;
	MPI_Offset FILESIZE;

        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
	// in SAIGE: FILESIZE = pvalVec.size();
	FILESIZE = 16; // number of rows
  	mpi_chunk_size = FILESIZE/size; //(FILESIZE-1)/size + 1; // number of rows assigned to each process (check for case when FILESIZE%size!=0)
  	mpi_k_start = rank*mpi_chunk_size;  // start row index for each process
  	//mpi_k_end =  std::min(FILESIZE, (rank+1)*mpi_chunk_size);
	//
	//allocate memory for the buffer holding the array (scatterd from the root process) local to each process
	int *local_chrArray = (int*)malloc(sizeof(int)*mpi_chunk_size);
	int *local_posArray = (int*)malloc(sizeof(int)*mpi_chunk_size);
	// send chunks of the vector to individual processes (to replicate the SAIGE formulation, although there the chunks are vectors and not arrays)
	if (rank == 0){
		// each of the vector below represent a column in the output file
		// these vectors are to be scattered to multiple processes	
		std::vector<int> chrVec{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
		std::vector<int> posVec{100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115};
		MPI_Scatter(chrVec.data(), mpi_chunk_size, MPI_INT, local_chrArray, mpi_chunk_size, MPI_INT, 0, MPI_COMM_WORLD);	
		MPI_Scatter(posVec.data(), mpi_chunk_size, MPI_INT, local_posArray, mpi_chunk_size, MPI_INT, 0, MPI_COMM_WORLD);	
	}
	else
	{
		MPI_Scatter(NULL, mpi_chunk_size, MPI_INT, local_chrArray, mpi_chunk_size, MPI_INT, 0, MPI_COMM_WORLD);	
		MPI_Scatter(NULL, mpi_chunk_size, MPI_INT, local_posArray, mpi_chunk_size, MPI_INT, 0, MPI_COMM_WORLD);	
	}
	// (SANITY CHECK) Get the average of the numbers in each process
	float local_avg;
	float local_sum=0;
	for(int i=0; i < mpi_chunk_size; i++)
		local_sum += local_chrArray[i];

	std::cout << "from process: " << rank << " of: " << size << " having elements " << mpi_chunk_size << " with local sum " << local_sum << std::endl;

	if (write_file_in_parallel == 1){	
	// write in parallel to a file
	// the 2 vectors are located at different locations in memory so they can't be written out directly in the correct format
	// One option is to create a 2D array from the different local vectors and then write out (creation of the 2D array for millions of rows can be expensive?)
	// The other option is to check how to write out non-contiguous arrays in parallel fashion
	
	// Option 1
	int nrows_global = 16; //total number of rows in the global vectors
	int ncols = 2;
	int chars_per_num = 9; // check this value
        MPI_Datatype num_as_string; //custom datatype where numbers are represented by strings
        MPI_Datatype localarray;    //custom datatype to store 2D arrays
        char *const fmt="%8.3f ";
        char *const endfmt="%8.3f\n";
	MPI_File file;


	float data[mpi_chunk_size][ncols];
        // fill the 2d array with the elements of the 2 1d arrays
        for(int i=0; i<ncols; i++){// loop over columns
            for(int j=0; j<mpi_chunk_size; j++){  // loop over local rows
                data[j][0] = local_chrArray[j]; // just use .at() when c++ vectors are being used
                data[j][1] = local_posArray[j];
            }
        }


        // creation of new MPI data type to represent numberes as ASCII chars
	// each number is represented by chars_per_num number of chars
	// MPI_Type_contagious defines a new data type that is a concatenation of a number of elements of an existing data type.
        // These replications are created into contiguous locations, resulting in a new contiguous data type.
        MPI_Type_contiguous(chars_per_num, MPI_CHAR, &num_as_string);
        MPI_Type_commit(&num_as_string);

        // to do: sprint data as char array and write using MPI-IO with Datatype MPI_CHAR

        // convert data into txt
	// sprintf: Composes a string with the same text that would be printed if format was used on printf, but instead of being printed, the content is stored as a C string in the buffer pointed by first argument (a string).
        char *data_as_txt = (char*)malloc(mpi_chunk_size*ncols*chars_per_num*sizeof(char));
        int count = 0; // keeping track of the start location in the buffer for storing the 2D array elements
        for (int i=0; i<mpi_chunk_size; i++) {
            for (int j=0; j<ncols-1; j++) {
                sprintf(&data_as_txt[count*chars_per_num], fmt, data[i][j]);
                count++;
            }
            sprintf(&data_as_txt[count*chars_per_num], endfmt, data[i][ncols-1]);
            count++;
        }

        std::cout << "rank: "<< rank << " data_as_txt: "<< data_as_txt << std::endl;

        // create a custom type for representing the local 2D arrays
        int global_sizes[2] = {nrows_global, ncols};
        int local_sizes [2] = {mpi_chunk_size, ncols};
        int starts[2]       = {mpi_k_start, 0}; //starting coordinates of the subarray in each dimension
        int order           = MPI_ORDER_C;

    	// MPI_Type_create_subarray creates an MPI datatype representing a subset of an array.
    	MPI_Type_create_subarray(2, global_sizes, local_sizes, starts, order, num_as_string, &localarray);
    	MPI_Type_commit(&localarray);
/*
	if (rank == 0){
        std::string filename("combined_data_header.txt");
        std::ofstream outfile;
        outfile.open(filename);
        outfile << "CHR\tPOS\t";
	outfile << "\n";
	outfile.close();
	}
*/
    	// open the file, and set the view
    	MPI_File_open(MPI_COMM_WORLD, "combined_data.txt",
                  MPI_MODE_CREATE|MPI_MODE_WRONLY,
                  //MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &file);

    	MPI_File_set_view(file, 0,  MPI_CHAR, localarray,
                           "native", MPI_INFO_NULL);

	std::cout << "LLL153" << std::endl;
    	MPI_File_write_all(file, data_as_txt, mpi_chunk_size*ncols, num_as_string, &status);
	std::cout << "LLL155" << std::endl;
    	MPI_File_close(&file);
	std::cout << "LLL157 from " << rank << std::endl;

    	MPI_Type_free(&localarray);
    	MPI_Type_free(&num_as_string);
	std::cout << "LLL161 from" << rank << std::endl;

    	//free(data[0]);
    	//free(data);
	}

	if (write_file_from_root == 1){
	// gather whole array at proc_0
	  if (rank == 0){
		int *gathered_chrArray = (int*)malloc(sizeof(int)*FILESIZE);
		int *gathered_posArray = (int*)malloc(sizeof(int)*FILESIZE);
		MPI_Gather(local_chrArray, mpi_chunk_size, MPI_INT, gathered_chrArray, mpi_chunk_size, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Gather(local_posArray, mpi_chunk_size, MPI_INT, gathered_posArray, mpi_chunk_size, MPI_INT, 0, MPI_COMM_WORLD);
		//std::cout << "gathered_vector " << gathered_vector[10] << std::endl;
		// write to a file
 		std::string filename("gathered_vector_proc_0.txt");
  		std::ofstream outfile;
  		outfile.open(filename);
		outfile << "CHR\tPOS\t";
		outfile << "\n";
		// loop over the vector elements
		for(int k=0; k<FILESIZE; k++){
			outfile << gathered_chrArray[k];
		  	outfile << "\t";
			outfile << gathered_posArray[k];
		  	outfile << "\n";
		}
		//outfile.close();

	  }
	  else{
		MPI_Gather(local_chrArray, mpi_chunk_size, MPI_INT, NULL, mpi_chunk_size, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Gather(local_posArray, mpi_chunk_size, MPI_INT, NULL, mpi_chunk_size, MPI_INT, 0, MPI_COMM_WORLD);
	  }
	}


	if (write_files_from_each_process == 1){
	// write separate file from each process
	  std::string filename("mpi_");
  	  filename += std::to_string(rank);
  	  filename += ".txt";

  	  std::ofstream mpi_outfile;
  	  mpi_outfile.open(filename);

	  // write out the header from proc_0
	  if (rank == 0)
		  //mpi_outfile << "CHR\tPOS\tMarkerID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\t";
		  mpi_outfile << "CHR\tPOS\t";
		  mpi_outfile << "\n";
	  //for(int k = mpi_k_start; k < mpi_k_end; k++) // loop over all the elements in the proc
	  for(int k = 0; k < mpi_chunk_size; k++){ // loop over all the elements in the proc
		  mpi_outfile << local_chrArray[k];
		  mpi_outfile << "\t";
		  mpi_outfile << local_posArray[k];
		  mpi_outfile << "\n";
	  }
	}
	std::cout << "LLL210 from" << rank << std::endl;
        free(local_chrArray);
        free(local_posArray);
	std::cout << "LLL213 from" << rank << std::endl;
	
        MPI_Finalize();
	std::cout << "LLL216 from" << rank << std::endl;

        return 0;

}


