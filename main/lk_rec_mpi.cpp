/* -*- -*- ----------------------------------------------------------
   TEPPP: Topological Entanglement in Polymers, Proteins and Periodic structures
   https://github.com/TEPPP-software/TEPPP.git
   Eleni Panagiotou, epanagio@asu.edu

   Copyright (2021) Eleni Panagiotou This software is distributed under
   the BSD 3-Clause License.

   See the README file in the top-level TEPPP directory.
   Contributors: Tom Herschberg, Kyle Pifer and Eleni Panagiotou
------------------------------------------------------------------------- */

#include "../include/funcs.h"
#include "mpi.h"

using namespace std;


/* -*- -*- ----------------------------------------------------------
   Takes as input the filename, the number of chains, the length of the chains, integer 0/1 statement (depending on whether the chains are closed, i.e. rings, or open, i.e. linear, 0: ring, 1: linear) and the box dimension (optional)
   If the box dimension is not specified, or if it is equal to 0, then the system is not periodic and the coordinates are unwrapped. 
   If a non-zero box-dimension is specified, the coordinates are unwrapped, according to the PBC.
   Returns the Gauss linking integral between each pair of chains
   Last line returns the mean absolute linking number over all pairs
------------------------------------------------------------------------- */


int main(int argc, char* argv[])
{
        MPI_Init(&argc, &argv);
        int rank, size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int num_species = 1;
	int total_num_chains;

	int num_chains = stoi(argv[5]);
	int chain_length = stoi(argv[6]);
	int ringlinear = stoi(argv[7]);

	total_num_chains = num_chains;

	// Arrays to store chain info for each species
	vector<int> chains_per_species;
	vector<int> lengths_per_species;
	vector<int> ringlinear_per_species;

	// Store first species info
	chains_per_species.push_back(num_chains);
	lengths_per_species.push_back(chain_length);
	ringlinear_per_species.push_back(ringlinear);

	if (argc > 9) {
		int num_chains_2 = stoi(argv[8]);
		int chain_length_2 = stoi(argv[9]);
		int ringlinear_2 = stoi(argv[10]);
		total_num_chains = total_num_chains + num_chains_2;
		num_species = num_species + 1;
		
		// Store second species info
		chains_per_species.push_back(num_chains_2);
		lengths_per_species.push_back(chain_length_2);
		ringlinear_per_species.push_back(ringlinear_2);
		cout<<"tototal 2 chains  "<<total_num_chains<<endl;

	}
	if (argc > 12) {
		int num_chains_3 = stoi(argv[11]);
		int chain_length_3 = stoi(argv[12]);
		int ringlinear_3 = stoi(argv[13]);
		total_num_chains = total_num_chains + num_chains_3;
		num_species = num_species + 1;
		
		// Store third species info
		chains_per_species.push_back(num_chains_3);
		lengths_per_species.push_back(chain_length_3);
		ringlinear_per_species.push_back(ringlinear_3);
	}
	if (argc > 15) {
		int num_chains_4 = stoi(argv[14]);
		int chain_length_4 = stoi(argv[15]);
		int ringlinear_4 = stoi(argv[16]);
		total_num_chains = total_num_chains + num_chains_4;
		num_species = num_species + 1;
		
		// Store fourth species info
		chains_per_species.push_back(num_chains_4);
		lengths_per_species.push_back(chain_length_4);
		ringlinear_per_species.push_back(ringlinear_4);
	}

	double box_dim_x;
	double box_dim_y;
	double box_dim_z;

    bool is_closed;
    
	box_dim_x = stod(argv[2]);
	box_dim_y = stod(argv[3]);
	box_dim_z = stod(argv[4]);

	cout<<box_dim_x<<" "<<box_dim_y<<" "<<box_dim_z<<endl;
	cout<<"tototal 0 chains  "<<total_num_chains<<endl;

	int num;
	double** coords;
	if (box_dim_x == 0)
	{
		coords = read_coords(argv[1], &num);
	}
	else
	{
		coords = read_coords_rec(argv[1], &num, chain_length, box_dim_x, box_dim_y, box_dim_z);
	}

    int chunk = total_num_chains / size;

        double result[chunk][total_num_chains];
        create_output_dir();
        string file_name = to_string(chain_length) + "_lk_rec_mpi_out_" + to_string(rank) + ".txt";
        ofstream outfile;
        if (!fs::exists("./output/lk_rec_mpi"))
        {
                cout << "Creating lk_rec_mpi directory..." << endl;
                fs::create_directory("./output/lk_rec_mpi");
        }
        outfile.open("./output/lk_rec_mpi/" + file_name);
	//ofstream outfile;
	int count = 0;
	double sum = 0;
	
	// Track is_closed for each species
	vector<bool> is_closed_per_species;
	for (int i = 0; i < num_species; i++) {
		is_closed_per_species.push_back(ringlinear_per_species[i] == 0);
	}
	
	// Calculate starting indices for each chain
	vector<int> chain_start_indices;
	int current_index = 0;
	for (int species = 0; species < num_species; species++) {
		for (int chain = 0; chain < chains_per_species[species]; chain++) {
			chain_start_indices.push_back(current_index);
			current_index += lengths_per_species[species];
		}
	}
	

	for (int chain_i = rank*chunk; chain_i < (rank+1)*chunk; chain_i++)
	{
		result[chain_i - (rank * chunk)][chain_i] = 0;
		// Determine which species this chain belongs to
		int species_i = 0;
		int chain_in_species_i = chain_i;
		for (int sp = 0; sp < num_species; sp++) {
			if (chain_in_species_i < chains_per_species[sp]) {
				species_i = sp;
				break;
			}
			chain_in_species_i -= chains_per_species[sp];
		}
		
		int length_i = lengths_per_species[species_i];
		bool is_closed_i = is_closed_per_species[species_i];
		
		// Extract chain_i coordinates
		double** chain1 = new double*[length_i];
		for (int j = 0; j < length_i; j++) {
			chain1[j] = new double[3];
			// Calculate the absolute index in the trajectory for this atom
			int atom_index = chain_start_indices[chain_i] + j;
			
			chain1[j][0] = coords[atom_index][0];
			chain1[j][1] = coords[atom_index][1];
			chain1[j][2] = coords[atom_index][2];
			
			// Optionally, print atom indices for debugging
			// cout << "Chain " << chain_i << ", atom " << j << " is at index " << atom_index << endl;
		}

		/**
		 * @brief This for loop is used to calculate the coordinates
		 * for both chain 1 and chain 2 for the lk calculation.
		 *
		 */
		// Loop through all other chains to form pairs
		for (int chain_j = chain_i + 1; chain_j < total_num_chains; chain_j++) {
			// Determine which species the second chain belongs to
			int species_j = 0;
			int chain_in_species_j = chain_j;
			for (int sp = 0; sp < num_species; sp++) {
				if (chain_in_species_j < chains_per_species[sp]) {
					species_j = sp;
					break;
				}
				chain_in_species_j -= chains_per_species[sp];
			}
			
			int length_j = lengths_per_species[species_j];
			bool is_closed_j = is_closed_per_species[species_j];
			
			// Extract chain_j coordinates
			double** chain2 = new double*[length_j];
			for (int k = 0; k < length_j; k++) {
				chain2[k] = new double[3];
				// Calculate the absolute index in the trajectory for this atom
				int atom_index = chain_start_indices[chain_j] + k;
				
				chain2[k][0] = coords[atom_index][0];
				chain2[k][1] = coords[atom_index][1];
				chain2[k][2] = coords[atom_index][2];
				
				// Optionally, print atom indices for debugging
				// cout << "Chain " << chain_j << ", atom " << k << " is at index " << atom_index << endl;
			}
			
			// Calculate linking number - make sure the lk function can handle different chain lengths
			// and different is_closed values for different chains
			double res = lk(chain1, chain2, length_i, length_j, is_closed_i && is_closed_j);
			outfile << "Chain " << chain_i << " and Chain " << chain_j << ": " << res << "\n";
			sum += abs(res);
			count++;

			result[chain_i - (rank * chunk)][chain_j] = res;
			// Clean up memory
			delete_array(chain2, length_j);
		}
		
		// Clean up memory
		delete_array(chain1, length_i);
	}

	delete_array(coords, num_chains);
	outfile.close();
	cout << "Absolute avg lk: " << sum / count << "\n";
        MPI_Finalize();
	return 0;
}
