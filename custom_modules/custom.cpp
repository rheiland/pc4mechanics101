/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

void setup_microenvironment( void )
{
	// set domain parameters
	
	initialize_microenvironment(); 	
	
	// these will ***overwrite*** values specified in the 
	// microenvironment_setup part of the XML,
	// based on what's in the user_parameters section 
	
	microenvironment.name = "synthetic tissue"; 
	
	// display the microenvironment again 
	
	microenvironment.display_information( std::cout ); 
	
	return; 
}

void create_cell_types( void )
{
	SeedRandom( parameters.ints("random_seed") ); 
	// housekeeping 
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL;   
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 

	cell_defaults.phenotype.mechanics.attachment_elastic_constant = parameters.doubles( "elastic_coefficient" );
	
	// static int cell_ID = find_cell_definition( "default" )->type; 
	
	Cell_Definition* pCD = find_cell_definition( "default" ); 
	pCD->functions.update_phenotype = pheno_update; 
	// pCD->functions.custom_cell_rule = NULL; // extra_elastic_attachment_mechanics;  // pre-1.8.0
	pCD->functions.custom_cell_rule = custom_cell_update;   // dt_mechanics
	pCD->functions.contact_function = standard_elastic_contact_function; 
	// pCD->functions.update_migration_bias = worker_cell_motility;
	pCD->phenotype.mechanics.attachment_elastic_constant = parameters.doubles( "elastic_const" );

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void pheno_update( Cell* pCell , Phenotype& phenotype , double dt )
{
	// return; 
	std::vector<Cell*> nearby = pCell->cells_in_my_container(); 
	
	// if at least 2 neighbors, turn off secretion 
		// if size >= 3, then we have "self" and at least two more 
	if( nearby.size() > 2 )
	{
		pCell->phenotype.secretion.set_all_secretion_to_zero(); 
		pCell->custom_data[ "secreting" ] = 0.0; 
		
		pCell->functions.update_phenotype = NULL; 
	}
	
	return; 
}

// BM adhesion-repulsion model (every dt_mechanics)
void custom_cell_update( Cell* pCell , Phenotype& phenotype , double dt )
{
    static int num_attached_BM = 0;
    static int attach_lifetime_i = pCell->custom_data.find_variable_index( "attach_lifetime" ); 
    static int attach_time_i = pCell->custom_data.find_variable_index( "attach_time" ); 
    static int attach_to_BM_i = pCell->custom_data.find_variable_index( "attach_to_BM" ); 

	// std::vector<Cell*> nearby = pCell->cells_in_my_container(); 
	
    // cap letters (X,N, etc) represent vectors
    // Agent vars:
    //   position = X
    //   radius = r
    //   adhesion radius = r_A

    if (pCell->ID == 0)
    {
        std::cout << "ID=0, phenotype.geometry.radius = " << phenotype.geometry.radius << std::endl;
    }

    double adhesion_radius = phenotype.geometry.radius * phenotype.mechanics.relative_maximum_adhesion_distance;
    int ncells_attached = 0;

    double displacement = 0.0 - pCell->position[1];  // displacement: just (negative) y (height) for test case
    //===================================
    //   attach
    //===================================
	if( pCell->custom_data[attach_to_BM_i] == 0.0 )  // not attached to BM
	{
        // double d = 0.0 - pCell->position[1];  // just (negative) y (height) for test case
        // std::cout << "t="<<PhysiCell_globals.current_time << ", ID=" << pCell->ID << ": d= " << d <<", adhesion radius= " << adhesion_radius << std::endl;
//        double pv = <0,-1,0> 
//        double nv = pv - d*nv;
        if (displacement <= 0.0 && displacement > -adhesion_radius )
        {
            std::cout << "t="<<PhysiCell_globals.current_time << "attaching ID=" << pCell->ID << ": displacement= " << displacement <<", adhesion radius= " << adhesion_radius << std::endl;
            // double p_BM = pv - d*nv
            pCell->custom_data[attach_to_BM_i] = 1.0;   // attached to BM now
            pCell->custom_data[attach_time_i] = 0.0;   // reset its time of being attached
            ncells_attached++;
            num_attached_BM++;
            // phenotype.motility.is_motile = false; 
        }
	}
    // else
    // {
    //     ncells_attached++;
    // }
    // displacement = X_BM - X = D

    pCell->custom_data[attach_time_i] += dt;

    //===================================
    //   mechanics
    //===================================
    pCell->velocity[1] += pCell->custom_data[attach_lifetime_i] * displacement;
    // axpy( &(pActingOn->velocity) , pao.mechanics.attachment_elastic_constant , displacement );   // cancer_immune_3D
	
    //===================================
    //   detach
    //===================================
    // if( UniformRandom() < dt / ( pCell->custom_data[attach_lifetime_i] + 1e-15 ) )
    if( pCell->custom_data[attach_time_i] > pCell->custom_data[attach_lifetime_i] )
    {
        pCell->custom_data[attach_to_BM_i] = 0.0;   // detach from BM 
        // detach_cells( pCell, pCell->state.attached_cells[0] ); 
        // phenotype.motility.is_motile = true; 
    }

    // if (pCell->ID == 51)
    //     std::cout << "(at ID=51)----- num_attached_BM = "<<num_attached_BM << std::endl;
    // std::cout << "---- custom_cell_update(): ncells_attached = " << ncells_attached << std::endl;
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	std::string color = "black"; 
	std::vector< std::string > output( 4 , color ); 
	
	// black cells if necrotic 
	if( pCell->phenotype.death.dead == true )
	{ return output; }

	output[3] = "none"; // no nuclear outline color 
	
	// static std::string worker_color = parameters.strings( "worker_color" ); 
	// static std::string cargo_color = parameters.strings( "cargo_color" ); 
	// static std::string director_color = parameters.strings( "director_color" ); 
	
	// static int worker_ID = find_cell_definition( "worker cell" )->type; 
	// static int cargo_ID = find_cell_definition( "cargo cell" )->type; 
	// static int director_ID = find_cell_definition( "director cell" )->type; 

	// if( pCell->type == worker_ID )
	// { color = worker_color; }
	// else if( pCell->type == cargo_ID )
	// { color = cargo_color; }
	// else if( pCell->type == director_ID )
	// { color = director_color; }
	// else
	// { color = "white"; } 
	
	// output[0] = color; 
	// output[2] = color; 
	
	return output; 
}

void setup_tissue( void )
{
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 		
	
	PhysiCell_SVG_options.length_bar = 200; 
	SVG_plot( "initial.svg" , microenvironment, 0.0 , 0.0 , my_coloring_function );	
	
	return; 
}

/*
void attach_cells( Cell* pCell_1, Cell* pCell_2 )
{
	#pragma omp critical
	{
		
	bool already_attached = false; 
	for( int i=0 ; i < pCell_1->state.neighbors.size() ; i++ )
	{
		if( pCell_1->state.neighbors[i] == pCell_2 )
		{ already_attached = true; }
	}
	if( already_attached == false )
	{ pCell_1->state.neighbors.push_back( pCell_2 ); }
	
	already_attached = false; 
	for( int i=0 ; i < pCell_2->state.neighbors.size() ; i++ )
	{
		if( pCell_2->state.neighbors[i] == pCell_1 )
		{ already_attached = true; }
	}
	if( already_attached == false )
	{ pCell_2->state.neighbors.push_back( pCell_1 ); }

	}

	return; 
}
*/

/*
void dettach_cells( Cell* pCell_1 , Cell* pCell_2 )
{
	#pragma omp critical
	{
		bool found = false; 
		int i = 0; 
		while( !found && i < pCell_1->state.neighbors.size() )
		{
			// if cell 2 is in cell 1's list, remove it
			if( pCell_1->state.neighbors[i] == pCell_2 )
			{
				int n = pCell_1->state.neighbors.size(); 
				// copy last entry to current position 
				pCell_1->state.neighbors[i] = pCell_1->state.neighbors[n-1]; 
				// shrink by one 
				pCell_1->state.neighbors.pop_back(); 
				found = true; 
			}
			i++; 
		}
	
		found = false; 
		i = 0; 
		while( !found && i < pCell_2->state.neighbors.size() )
		{
			// if cell 1 is in cell 2's list, remove it
			if( pCell_2->state.neighbors[i] == pCell_1 )
			{
				int n = pCell_2->state.neighbors.size(); 
				// copy last entry to current position 
				pCell_2->state.neighbors[i] = pCell_2->state.neighbors[n-1]; 
				// shrink by one 
				pCell_2->state.neighbors.pop_back(); 
				found = true; 
			}
			i++; 
		}

	}
	
	return; 
}
*/

/*
void add_elastic_velocity( Cell* pActingOn, Cell* pAttachedTo , double elastic_constant )
{
	std::vector<double> displacement = pAttachedTo->position - pActingOn->position; 
	axpy( &(pActingOn->velocity) , elastic_constant , displacement ); 
	
	return; 
}
*/

/*
void extra_elastic_attachment_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	// if I am 
	std::vector<double> velocity(3,0.0); 
	
	for( int i=0; i < pCell->state.neighbors.size() ; i++ )
	{
		add_elastic_velocity( pCell, pCell->state.neighbors[i], pCell->custom_data["elastic coefficient"] ); 
	}

	return; 
}	
*/
