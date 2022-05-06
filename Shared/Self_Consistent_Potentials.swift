//
//  Self_Consistent_Potentials.swift
//  SelfConsistent_Hartree_Fock_Slater_non_relativistic
//
//  Created by Varrick Suezaki on 8/11/18.
//  Copyright © 2018 Varrick Suezaki. All rights reserved.
//

import Foundation

class Self_Consistent_Potentials: NSObject {
    
    var myHartreeFockSCFCalculator: HermanSkillmanCalculator? = nil
    
    /// Name: calculate_total_charge_density
    /// Description: Calculates total normalized charge density from number of electrons per orbital and P wave function values
    /// Equation: σ(r) = -Σω|P(r)|^2 (ω is occupation number)
    ///
    /// - Parameters:
    ///   - total_mesh_count: Total number of points in mesh
    ///   - core_shells: Number of core shells
    ///   - val_shells: Number of valence shells
    ///   - results_array: Array of containing r, P wave function, quantum numbers, number of electrons and energy values for each orbital
    /// - Returns: Returns total normalized charge density of current orbital as a function of r
    func calculate_total_charge_density(total_mesh_count: Int, core_shells: Int, val_shells: Int, results_array: [(r_list: [Double], psi_list: [Double], quant_n: Double, quant_l: Double, quant_m: Double, number_electrons: Double, new_energy: Double)]) -> [Double] {
        
        // Initializes normalized charge density arrays for core electrons and valence electrons
        var rcores_array: [Double] = Array(repeating: 0.0, count: total_mesh_count)
        var rvals_array: [Double] = Array(repeating: 0.0, count: total_mesh_count)
        
        // Array for total normalized charged density values
        var rtot_array: [Double] = []
        
        // Total number of shells
        let total_shells: Int = core_shells + val_shells
        
        // For loop that loops over total number of shells.  While the count is not larger than the number of core shells, the normalized charged density is calculated for the core shells, else its calculated for the valence shells
        for i in 1...total_shells {
            
            if i <= core_shells {
                
                for j in 0..<results_array[i-1].psi_list.count {
                    rcores_array[j] = rcores_array[j] + results_array[i-1].number_electrons*pow(results_array[i-1].psi_list[j], 2.0)
                }
                
            } else {
                
                for j in 0..<results_array[i-1].psi_list.count {
                    rvals_array[j] = rvals_array[j] + results_array[i-1].number_electrons*pow(results_array[i-1].psi_list[j], 2.0)
                }
                
            }
            
        }
        
        // Adds each element of core and valence charge density arrays for total charge density array
        for i in 0..<total_mesh_count {
            rtot_array.append(rcores_array[i] + rvals_array[i])
        }
        
        return(rtot_array)
        
    }
    
    /// Name: modified_potential
    /// Description: Calculates modified Hartree-Fock-Slater potential
    /// Equation: V0(r) = -2*Z/r - (2/r)⎰(σ(t)*dt) - 2*⎰((σ(t)/t)*dt) - 6*[(3/8π)*|ρ(r)|]^(1/3)
    ///
    /// - Parameters:
    ///   - pot_r_list: List of r values 
    ///   - rtot_array: Total normalized charge density array calculated by calculate_total_charge_density function
    ///   - z_value: Z value of atom
    ///   - exchange_alpha: Scalar to multiply exchange term by.
    ///   - ionticity: Net charge of atom
    /// - Returns: Returns modified Hartree-Fock-Slater potential array
    func modified_potential(pot_r_list: [Double], rtot_array: [Double], z_value: Double, exchange_alpha: Double, ionticity: Double) -> [Double] {

        // Instance of Functional_functions class
        let potential_functional_functions_inst = myHartreeFockSCFCalculator!.functional_functions_inst
        
        // Initialize first integral, second integral and final modified potential array
        var second_rtot_array: [Double] = []
        var pot_final: [Double] = []
        let first_rtot_array: [Double] = rtot_array.map({$0})
        
        // Elements for second integral
        for i in 0..<pot_r_list.count {
            second_rtot_array.append(rtot_array[i]/pot_r_list[i])
        }
        
        // Loops over each r value and calculates first integral, second integral and exchange term elements
        for i in 0..<pot_r_list.count-1 {
            // exchange term
            let pot_exchange_value: Double = -6.0*exchange_alpha*pow((3.0*pot_r_list[i]*first_rtot_array[i])/315.82734, 1.0/3.0)
            
            // first integral
            let first_integral: Double = 2.0*potential_functional_functions_inst.integrate_mesh(x_list: pot_r_list, y_list: first_rtot_array, lower_bound_index: 0, upper_bound_index: i)

            // second integral
            let second_integral: Double = 2.0*pot_r_list[i]*potential_functional_functions_inst.integrate_mesh(x_list: pot_r_list, y_list: second_rtot_array, lower_bound_index: i, upper_bound_index: pot_r_list.count-1)

            // Adds them all together and adds it to pot_final array
            pot_final.append(-2.0*z_value + first_integral + second_integral + pot_exchange_value)
        }
        
        // If any of the final potential values are larger than -2, automatically set the values to -2.  The point where the final potential values become larger than -2 is an indicator to correct the potential to the right aymptotic behavior
        for i in 0..<pot_final.count {
            
            if pot_final[i] > -2.0 {
                
                pot_final[i] = -2.0
                
            }
            
        }
        
        // All terms besides nuclear term should be 0
        pot_final[0] = -2.0*z_value
        
        return(pot_final)
        
    }
    /// Name: interpolate_Potential_And_r
    /// Description: Takes r_list and potential and converts them to a consistent point density, which is that of the start of r_list
    func interpolate_Potential_And_r(initial_r_list: [Double], initial_Potential: [Double])->(r_list: [Double], Potential: [Double], codex: [(initial: Int,interpolated: Int)]){
        let r_step = initial_r_list[1] - initial_r_list[0]
        let r_max = initial_r_list.max()!
        let number_of_steps = round(r_max/r_step)
        var index_codex_for_de_interpolation: [(initial: Int,interpolated: Int)] = []

        print("r_step = \(r_step), r_max = \(r_max), number_of_steps = \(number_of_steps)")
        var new_r_list: [Double] = []
        var new_potential_list: [Double] = []

        for i in stride(from: 0, through: number_of_steps, by: 1){
            new_r_list.append(Double(i)*r_step)
        }
        
        new_potential_list.append(initial_Potential[0])
        
        print("initial_r_list.count - 1 = \(initial_r_list.count - 1 )")
        for i in stride(from: 0, to: initial_r_list.count-1, by: 1){
            let r_dif = initial_r_list[i+1] - initial_r_list[i]
            let steps_between = r_dif/r_step
            let v_step = (initial_Potential[i+1]-initial_Potential[i])/steps_between
            let v0 = initial_Potential[i]
//            if(Int(steps_between) > 1){
                
                for j in stride(from: 1, through: Int(round(steps_between)), by: 1){
                    if(j == 1){
                        index_codex_for_de_interpolation.append((initial: i, interpolated: new_potential_list.count-1))
                    }
                    new_potential_list.append(Double(j)*v_step+v0)
                    
                }
//            }
        
        }

//        if(new_potential_list.count < new_r_list.count){
//            let number_to_append = new_r_list.count - new_potential_list.count
//            for i in stride(from: 1, through: number_to_append, by: 1){
//                new_potential_list.append(initial_Potential.last!)
//            }
//        }
        print("new_r_list is \(new_r_list.count) long and new_potential_list is \(new_potential_list.count) long")
        return((r_list: new_r_list, Potential: new_potential_list, codex: index_codex_for_de_interpolation))
    }
    
    /// Name: initial_input_potential
    /// Description: Takes input 110 potential points and expands it to 441
    ///
    /// - Parameters:
    ///   - KEY: Int that specifies input potential format KEY=0 only one implemented
    ///   - z_value: Z value of atom
    ///   - initial_input_pot_list: Initial input potential from input file
    ///   - ionticity: Net charge of atom
    ///   - mesh_count: Total number of points in mesh
    func initial_input_potential(KEY: Int, z_value: Double, initial_input_pot_list: [Double], ionticity: Double, mesh_count: Double) -> [Double]{
        
        var new_input_pot_list: [Double] = []
        let pot_last_index: Int = initial_input_pot_list.count-1
        
        // Adds two more values to initial_input_pot_list so when it's expanded, there are 441 points
        let modi_initial_input_pot_list: [Double] = initial_input_pot_list + [initial_input_pot_list[pot_last_index],initial_input_pot_list[pot_last_index]]
        
        var block_counter: Int = 9
        
        // Indicates what format input potential is in.  0 for 110 point, 1 for 441 point and 2 for 2 sets of 441 point potentials (This code assumes 0)
        switch KEY {
            
        case 0:
            
            // loops over each value in the 110 point input potential
            for i in 0..<pot_last_index+1 {
                
                block_counter -= 1
                
                // Boolean to indicate whether or not loop is on first block (expansion differs in beginning compared to exceeding blocks)
                let block_counter_bool: Bool = block_counter >= 0
                
                switch block_counter_bool {
                    
                    // Expands current point into 4 different points for all blocks except starting black
                case true:
                    
                    let r1_approx: Double = (21.0*modi_initial_input_pot_list[i] + 14.0*modi_initial_input_pot_list[i+1] - 3.0*modi_initial_input_pot_list[i+2])/32.0
                    let r2_approx: Double = (3.0*modi_initial_input_pot_list[i] + 6.0*modi_initial_input_pot_list[i+1] - modi_initial_input_pot_list[i+2])/8.0
                    let r3_approx: Double = (5.0*modi_initial_input_pot_list[i] + 30.0*modi_initial_input_pot_list[i+1] - 3.0*modi_initial_input_pot_list[i+2])/32.0
                    
                    new_input_pot_list.append(contentsOf: [modi_initial_input_pot_list[i],r1_approx,r2_approx,r3_approx])
                    
                    // Expands current point into 4 different points for beginning block only
                default:
                    
                    let last_r1_approx: Double = (22.0*modi_initial_input_pot_list[i] + 11.0*modi_initial_input_pot_list[i+1] - modi_initial_input_pot_list[i+2])/32.0
                    let last_r2_approx: Double = (10.0*modi_initial_input_pot_list[i] + 15.0*modi_initial_input_pot_list[i+1] - modi_initial_input_pot_list[i+2])/24.0
                    let last_r3_approx: Double = (6.0*modi_initial_input_pot_list[i] + 27.0*modi_initial_input_pot_list[i+1] - modi_initial_input_pot_list[i+2])/32.0
                    
                    new_input_pot_list.append(contentsOf: [modi_initial_input_pot_list[i],last_r1_approx,last_r2_approx,last_r3_approx])
                    
                    block_counter = 9
                    
                }
                
            }
            
            new_input_pot_list.append(initial_input_pot_list[pot_last_index])
            
            // Adds values to new potential array if there are more than 441 points
            if mesh_count > 441.0 {
                for _ in 442...Int(mesh_count) {
                    new_input_pot_list.append((ionticity + 1.0)/z_value)
                }
            }
            // Cases 1 and 2 are ignored since this code assumes input potential is in the 110 point form
        case 1:
            
            let _: Int = 0
            
        case 2:
            
            let _: Int = 0
            
        default:
            
            print("KEY input not 0, 1 or 2")
            
        }
        // Get the unnormalized potential
        new_input_pot_list = new_input_pot_list.map({-2.0*z_value*$0})
        print(new_input_pot_list)
        return(new_input_pot_list)
        
    }
    
    /// Name: unpack_r_and_wavefunction
    /// Description: reverts interpolated mesh to coarser mesh of ground ground states
    ///
    /// - Parameters:
    ///   - codex: tuple containing the corresponding indeces of each array
    ///   - r_values: fine r values t ounpack
    ///   - psi_values: psi values to unpack
    /// - Returns: A tuple containing the unpacked r and psi
    func unpack_r_and_wavefunction(codex: [(initial: Int,interpolated: Int)], r_values: [Double], psi_values: [Double]) -> (r_values: [Double], psi_values: [Double]){
        
        var unpacked_r: [Double] = []
        var unpacked_psi: [Double] = []
        
        for i in codex{
            unpacked_r.append(r_values[i.interpolated])
            unpacked_psi.append(psi_values[i.interpolated])
            
        }
        if (unpacked_r.count == 0 && unpacked_psi.count == 0){
            print("unpacked to zero")
            print(r_values.count, psi_values.count)
        }
        return((r_values: unpacked_r, psi_values: unpacked_psi))
    }
    
    
    /// Name: calculate_ionic_potential
    /// Description: Calculates modified input potential if it's an ion
    ///
    /// - Parameters:
    ///   - pot_r_list: List of r values
    ///   - pot_list: List of potential values
    ///   - ionticity: Net charge of atom
    ///   - branching_ratio: Ratio of the rate constant for a particular product
    ///   - ionic_radius: Radius of ion
    /// - Returns: Array of modified input potential
    func calculate_ionic_potential(pot_r_list: [Double], pot_list: [Double], ionticity: Double, branching_ratio: Double, ionic_radius: Double) -> [Double] {
        
        // Calculates ionticity taking into account the branching ratio
        let new_ionticity: Double = ionticity*branching_ratio
        var new_ionic_pot_list: [Double] = []
        var ionic_radius_bool: Bool = true
        
        // Loops over each value in pot_list (should be same as pot_r_list)
        for i in 0..<pot_r_list.count {
            
            // Boolean indicating whether current radius is larger or smaller than ionic radius
            ionic_radius_bool = pot_r_list[i] < ionic_radius
            
            // Switch case to calculate potential depending on current radius value
            switch ionic_radius_bool {
                
            case true:
                
                new_ionic_pot_list.append(pot_list[i] + 2.0*(new_ionticity*pot_r_list[i]/ionic_radius))
                
            default:
                
                new_ionic_pot_list.append(pot_list[i] + 2.0*new_ionticity)
                
            }
            
        }
        
        return(new_ionic_pot_list)
        
    }
    
    /// Name: self_consistent_potential
    /// Description: Calculates self-consistent potential according to parameter 'beta' criteria.  Self-consistent energy and P wave function values are calculated using the 'schroedinger_subroutine' function.  Potentials are re-calculated using the P wave function values and compared to the input potential.  The maximum difference between the two potential arrays must be equal to or smaller than beta to be considered self consistent.
    ///
    /// - Parameters:
    ///   - z_value: Z value of atom
    ///   - delta_x_initial: Initial step in x values
    ///   - number_of_blocks: Number of blocks in mesh
    ///   - scalar: Scalar to convert x values to r values
    ///   - number_points: Number of points in mesh
    ///   - user_input_pot_list: Initial input potential selected by user
    ///   - core_shells: Number of core electron shells
    ///   - val_shells: Number of valence electron shells
    ///   - electron_config_array: Array of initial input energies and quantum numbers for each orbital
    ///   - beta_criterion: Criterion that measures potential self-consistency
    ///   - exchange_alpha: Scalar to multiply exchange term by.
    ///   - pratt_alpha: Double ranging from 0.0-1.0 that determines amount potential are changed between iterations (0.7-0.9 is slow but is more likely to converge, 0.1-0.5 is fast but not likely to converge correctly).
    ///   - KEY: Indicates input potential format
    ///   - thresh_criterion: Acceptable value of deltaE/E
    ///   - max_beta_iterations: Max number of iterations to calculate self-consistent potentials
    ///   - max_thresh_iterations: Max number of iterations to calculate self-consistent energy and P wave function values
    ///   - ionic_radius: Radius of ion if input atom is ion
    ///   - branching_ratio: Branching ratio for ionic radius if atom is ion
    ///   - ionticity: Net charge of atom
    /// - Returns: Returns a tuple with three values, the first is an tuple array containing r, P wave function, energy, # elecrons and quantum numbers of each orbital.  The second is an array of the final self-consistent potential.  The third value is an array of the full r_mesh.
    func self_consistent_potential(z_value: Double, delta_x_initial: Double, number_of_blocks: Double, scalar: Double, number_points: Double, user_input_pot_list: [Double], core_shells: Int, val_shells: Int, electron_config_array: [(quant_n: Double, quant_l: Double, quant_m: Double, numb_electrons: Double, trial_energy: Double)], beta_criterion: Double, exchange_alpha: Double, pratt_alpha: Double, KEY: Int, thresh_criterion: Double, max_beta_iterations: Int, max_thresh_iterations: Int, ionic_radius: Double, branching_ratio: Double, ionticity: Double, photoionization_energies: [Double]) -> ([(r_list: [Double], psi_list: [Double], quant_n: Double, quant_l: Double, quant_m: Double, number_electrons: Double, new_energy: Double)], [Double], [Double]){
        
        // Instances of classes initialized in View Controller
        let potential_schrod_eq_subroutine_inst = myHartreeFockSCFCalculator!.schrod_eq_subroutine_inst
        let Photoionization_schrod_eq_subroutine_inst = myHartreeFockSCFCalculator!.schrod_eq_subroutine_inst
        let potential_wavefunction_values_inst = myHartreeFockSCFCalculator!.wavefunction_values_inst
        let potential_mesh_potential_init_inst = myHartreeFockSCFCalculator!.mesh_potential_init_inst
        
        print("pratt_alpha = \(pratt_alpha)")
        potential_schrod_eq_subroutine_inst.myHartreeFockSCFCalculator = myHartreeFockSCFCalculator
        
        // Initializes total mesh points and initial input potential array
        let total_mesh_count: Int = Int(number_points*number_of_blocks) + 1
        var input_pot_list: [Double] = self.initial_input_potential(KEY: KEY, z_value: z_value, initial_input_pot_list: user_input_pot_list, ionticity: ionticity, mesh_count: Double(total_mesh_count))
//        print("first input_pot_list = \(input_pot_list)")
        
        // Initializes variables to hold calculated beta values and beta loop counter.  Also initializes arrays to hold calculated values (mesh values, P wavefunction values etc.)
        var current_beta_max: Double = 1.0
        var beta_counter: Int = 0
        var results_array: [(r_list: [Double], psi_list: [Double], quant_n: Double, quant_l: Double, quant_m: Double, number_electrons: Double, new_energy: Double)] = []
        
        var r_list: [Double] = []
        var psi_list: [Double] = []
        var full_mesh: [Double] = []
        
        // Sets a range for the expected energy
        var left_energy_scalar: Double = 1.65
        var right_energy_scalar: Double = 0.35
        
        // If its an ion, change lower bound to increase range
        if ionticity > 0.0 {
            left_energy_scalar = 2.4
            right_energy_scalar = 0.35
        }
        
        // While loop that compares max calculated beta value to accepted beta criterion.  Modifies initial potential, and starts calculation over until criterion is met or max number of iterations is exceeded
        while (current_beta_max >= beta_criterion) || (current_beta_max == 0.0) {
//        }
            
//        while (abs(current_beta_max - beta_criterion) >= 1e-3 ) {
            
            results_array = []
            
            // Loops over each orbital
            for i in electron_config_array {
//                print("trial_energy = \(i.trial_energy)")
                r_list = []
                psi_list = []
                
                // Calculates P wavefunction, potential and orbital energy values for the ground state
//                print("number of blocks = \(number_of_blocks), number_points = \(number_points)")
//                print("input_pot_list for l = \(input_pot_list.count) long")
                potential_schrod_eq_subroutine_inst.schroedinger_subroutine(z_value: z_value, trial_energy: i.trial_energy, number_blocks: number_of_blocks, l_number: i.quant_l, number_points: number_points, initial_delta_x: delta_x_initial, mesh_scalar: scalar, principal_quant_number: i.quant_n, input_pot_list: input_pot_list, thresh_criterion: thresh_criterion, max_thresh_iterations: max_thresh_iterations, left_energy_scalar: left_energy_scalar, right_energy_scalar: right_energy_scalar)
                
                // Appends calculated radius and wavefunction values to r and psi list
                for i in potential_wavefunction_values_inst.norm_Pwavefunction_tuple_array {
                    for j in i.wavefunction_list {
                        r_list.append(j.r_value)
                        psi_list.append(j.wavefunction_value)
                    }
                }
                // Adds values to results array
//                print("l r_list is \(r_list.count) long and psi_list is \(psi_list.count) long")
                results_array.append((r_list: r_list, psi_list: psi_list, quant_n: i.quant_n, quant_l: i.quant_l, quant_m: i.quant_m, number_electrons: i.numb_electrons, new_energy: potential_wavefunction_values_inst.energy_value))
                
                full_mesh = potential_mesh_potential_init_inst.r_mesh
//                print(full_mesh)
                //Calculate Excited Wavefunctions for l+1 and l-1 for all n, m, and photon energy and append those results to their own structures
                //This takes a lot longer than before simply because each ground state has two additional wavefunctions for each energy, and there can be 20+ energies, or 40+ wavefunctions.
                //This would be a good thing to thread if possible.
                
                
                
                // Clears all values to start next orbital calculation
                potential_mesh_potential_init_inst.KE_term_list.removeAll()
                potential_mesh_potential_init_inst.r_mesh = [0.0]
                potential_mesh_potential_init_inst.trial_pot_list.removeAll()
                potential_mesh_potential_init_inst.delta_r_list.removeAll()
                potential_schrod_eq_subroutine_inst.test_final_values_list.removeAll()
                potential_schrod_eq_subroutine_inst.test_final_r_values_list.removeAll()
                potential_wavefunction_values_inst.outward_Pwavefunction_tuple_array.removeAll()
                potential_wavefunction_values_inst.inward_Pwavefunction_tuple_array.removeAll()
                potential_wavefunction_values_inst.inward_log_deriv_integral_values.removeAll()
                potential_wavefunction_values_inst.outward_log_deriv_integral_values.removeAll()
                potential_wavefunction_values_inst.norm_Pwavefunction_tuple_array.removeAll()
                
            }
            
            // Calculates total electron density using wavefunction values and orbital values.  Using electron density, Z value, exchange alpha constant, potential values and radius values, modified potential is calculated
            let rtot_array: [Double] = self.calculate_total_charge_density(total_mesh_count: total_mesh_count, core_shells: core_shells, val_shells: val_shells, results_array: results_array)
            let pot_r_list: [Double] = full_mesh
            var pot_final: [Double] = self.modified_potential(pot_r_list: pot_r_list, rtot_array: rtot_array, z_value: z_value, exchange_alpha: exchange_alpha, ionticity: ionticity)
//            print("pot_final = \(pot_final)")
            // If its an ion, uses a slightly different program to calculate final modified potential
            if ionic_radius > 0.0 {
                pot_final = self.calculate_ionic_potential(pot_r_list: pot_r_list, pot_list: pot_final, ionticity: ionticity, branching_ratio: branching_ratio, ionic_radius: ionic_radius)
            }
            
            // Initializes and calculates beta values.  Finds maximum value from array
            var beta_list: [Double] = [0.0]
            for i in 0..<pot_r_list.count-1 {
                beta_list.append(abs(input_pot_list[i] - pot_final[i]))
            }
            current_beta_max = beta_list.max()!
            
            // Compares max calculated beta value to beta criterion.  If it's larger, modify initial potential using pratt scheme (in herman skillman book).
            if (current_beta_max > beta_criterion) || (current_beta_max == 0.0) {
                
                for i in 0..<pot_r_list.count-1 {
                    
                    input_pot_list[i] = pratt_alpha*input_pot_list[i] + (1.0 - pratt_alpha)*pot_final[i]
                }
                
                print("CURRENT BETA: \(current_beta_max)")
                
                if current_beta_max == 0.0 {
//                    print("input pot list = ")
//                    print(input_pot_list)
                    for i in results_array {
                        print(i.new_energy)
                    }
                }
                
            }
            
            beta_counter += 1

    }
        
        print("")
        print("ITERATIONS: \(beta_counter) FINAL BETA: \(current_beta_max)")
        for i in results_array {
            print(i.new_energy)
        }
        print("")
        
        return(results_array, input_pot_list, full_mesh)
        
    }
    
    
    
    
    
    
    func self_consistent_potential2(z_value: Double, delta_x_initial: Double, number_of_blocks: Double, scalar: Double, number_points: Double, user_input_pot_list: [Double], core_shells: Int, val_shells: Int, electron_config_array: [(quant_n: Double, quant_l: Double, quant_m: Double, numb_electrons: Double, trial_energy: Double)], beta_criterion: Double, exchange_alpha: Double, pratt_alpha: Double, KEY: Int, thresh_criterion: Double, max_beta_iterations: Int, max_thresh_iterations: Int, ionic_radius: Double, branching_ratio: Double, ionticity: Double, photoionization_energies: [Double], r_mesh: [Double]) -> ([(r_list: [Double], psi_list: [Double], quant_n: Double, quant_l: Double, quant_m: Double, number_electrons: Double, photo_energy: Double, new_energy: Double)], [(r_list: [Double], psi_list: [Double], quant_n: Double, quant_l: Double, quant_m: Double, number_electrons: Double, photo_energy: Double, new_energy: Double)], [Double], [Double]){
            
        
        print("started second subroutine")
            // Instances of classes initialized in View Controller
            let potential_schrod_eq_subroutine_inst = myHartreeFockSCFCalculator!.schrod_eq_subroutine_inst
            let potential_wavefunction_values_inst = myHartreeFockSCFCalculator!.wavefunction_values_inst
            let potential_mesh_potential_init_inst = myHartreeFockSCFCalculator!.mesh_potential_init_inst
            
            potential_schrod_eq_subroutine_inst.myHartreeFockSCFCalculator = myHartreeFockSCFCalculator
            
            // Initializes total mesh points and initial input potential array
            let total_mesh_count: Int = Int(number_points*number_of_blocks) + 1
//            var input_pot_list: [Double] = self.initial_input_potential(KEY: KEY, z_value: z_value, initial_input_pot_list: user_input_pot_list, ionticity: ionticity, mesh_count: Double(total_mesh_count))
        
        
            //Interpolates r and potential meshes to increase precision
            let interpolated_r_and_v_tuple = interpolate_Potential_And_r(initial_r_list: r_mesh, initial_Potential: user_input_pot_list)
            let interpolated_Potential = interpolated_r_and_v_tuple.Potential
//            print(interpolated_Potential)
            let interpolated_r_mesh = interpolated_r_and_v_tuple.r_list
            let interpolation_codex = interpolated_r_and_v_tuple.codex
            
            // Initializes variables to hold calculated beta values and beta loop counter.  Also initializes arrays to hold calculated values (mesh values, P wavefunction values etc.)
            var current_beta_max: Double = 1.0
            var beta_counter: Int = 0
            var results_array: [(r_list: [Double], psi_list: [Double], quant_n: Double, quant_l: Double, quant_m: Double, number_electrons: Double, new_energy: Double)] = []
            var results_array_For_l_Minus_One: [(r_list: [Double], psi_list: [Double], quant_n: Double, quant_l: Double, quant_m: Double, number_electrons: Double, photo_energy: Double, new_energy: Double)] = []
            var results_array_For_l_Plus_One: [(r_list: [Double], psi_list: [Double], quant_n: Double, quant_l: Double, quant_m: Double, number_electrons: Double, photo_energy: Double, new_energy: Double)] = []
            var r_list: [Double] = []
            let input_mesh = r_mesh
            var psi_list: [Double] = []

        var unpackedTuple: (r_values: [Double], psi_values: [Double]) = (r_values: [], psi_values: [])
       
            var full_mesh: [Double] = []
            
            // Sets a range for the expected energy
            var left_energy_scalar: Double = 1.65
            var right_energy_scalar: Double = 0.35
            
            // If its an ion, change lower bound to increase range
            if ionticity > 0.0 {
                left_energy_scalar = 2.4
                right_energy_scalar = 0.35
            }
            
            // While loop that compares max calculated beta value to excepted beta criterion.  Modifies initial potential, and starts calculation over until criterion is met or max number of iterations is exceeded
//            while (current_beta_max >= beta_criterion) || (current_beta_max == 0.0) {
                
            //while (abs(current_beta_max - beta_criterion) >= 1e-3 ) {
                
                results_array = []
                
                // Loops over each orbital
        
        print("target number of wavefunctions should be \(electron_config_array.count*photoionization_energies.count)")
//        print(electron_config_array)
        
        print("photoionization_energies = \(photoionization_energies)")
                for i in electron_config_array {
                    
                    r_list = []
                    psi_list = []
//                    print(i.trial_energy)
                    
                    for photon_energy in photoionization_energies{
                        
                        var state_energy = photon_energy + i.trial_energy
//                        print("hv = \(photon_energy), ground_energy = \(i.trial_energy), kinetic_energy = \(state_energy)")
                        //l+1
                        r_list = []
                        psi_list = []
                     
                        if(state_energy < 0){
//                            print("state_energy = \(state_energy), n = \(Int(i.quant_n)), l+1 = \(Int(i.quant_l+1.0)) triggered negative state energy")
//                            print("l+1 number of blocks = \(number_of_blocks), number_points = \(number_points)")

                            
                            for i in stride(from: 0, to: interpolated_r_mesh.count, by: 1){
                                
                                r_list.append(interpolated_r_mesh[i])
                                psi_list.append(0.0)
                            }
                            
                        }
                        else{
                            if(state_energy < 0.003){
                                state_energy = 0.0
                            }
//                            print("state_energy = \(state_energy), n = \(Int(i.quant_n)), l+1 = \(Int(i.quant_l+1.0)) triggered positive state energy")
//                            print("input_pot_list for l+1 = \(input_pot_list.count) long")
    //                        print("l+1number of blocks = \(number_of_blocks), number_points = \(number_points)")
                            potential_schrod_eq_subroutine_inst.schroedinger_photoionization_subroutine(z_value: z_value, trial_energy: state_energy, number_blocks: number_of_blocks, l_number: i.quant_l+1, number_points: number_points, initial_delta_x: delta_x_initial, mesh_scalar: scalar, principal_quant_number: i.quant_n, input_pot_list: interpolated_Potential, thresh_criterion: thresh_criterion, max_thresh_iterations: max_thresh_iterations, left_energy_scalar: left_energy_scalar, right_energy_scalar: right_energy_scalar, r_list: interpolated_r_mesh)
                            
                            for i in potential_wavefunction_values_inst.norm_Pwavefunction_tuple_array {
                                for j in i.wavefunction_list {
                                    r_list.append(j.r_value)
                                    psi_list.append(j.wavefunction_value)
                                }
                            }
                        }
                        
                        //Unpack fine wavefunction to coarser grind for integration
                        unpackedTuple = unpack_r_and_wavefunction(codex: interpolation_codex, r_values: r_list, psi_values: psi_list)
                        if(unpackedTuple.r_values.count == 0 || unpackedTuple.psi_values.count == 0){
                            unpackedTuple.r_values.append(0.0)
                            unpackedTuple.psi_values.append(0.0)
                        }
                        // Adds values to results array
    //                    print("l+1 r_list is \(r_list.count) long and psi_list is \(psi_list.count) long")
//                        print("Photoion n = \(i.quant_n)")
                        results_array_For_l_Plus_One.append((r_list: unpackedTuple.r_values, psi_list: unpackedTuple.psi_values, quant_n: i.quant_n, quant_l: i.quant_l, quant_m: i.quant_m, number_electrons: i.numb_electrons, photo_energy: photon_energy, new_energy: state_energy-photon_energy))

                        full_mesh = unpackedTuple.r_values

                        potential_mesh_potential_init_inst.KE_term_list.removeAll()
                        potential_mesh_potential_init_inst.r_mesh = [0.0]
                        potential_mesh_potential_init_inst.trial_pot_list.removeAll()
                        potential_mesh_potential_init_inst.delta_r_list.removeAll()
                        potential_schrod_eq_subroutine_inst.test_final_values_list.removeAll()
                        potential_schrod_eq_subroutine_inst.test_final_r_values_list.removeAll()
                        potential_wavefunction_values_inst.outward_Pwavefunction_tuple_array.removeAll()
                        potential_wavefunction_values_inst.inward_Pwavefunction_tuple_array.removeAll()
                        potential_wavefunction_values_inst.inward_log_deriv_integral_values.removeAll()
                        potential_wavefunction_values_inst.outward_log_deriv_integral_values.removeAll()
                        potential_wavefunction_values_inst.norm_Pwavefunction_tuple_array.removeAll()
    //
                        //l-1
                        //l CANNOT BE NEGATIVE
                        r_list = []
                        psi_list = []
                        
                        if(i.quant_l-1.0 < 0.0){
                            r_list.append(0.0)
                            psi_list.append(0.0)
                        }
                        else{
                            
                            if(state_energy < 0){

                                for i in stride(from: 0, to: input_mesh.count, by: 1){
                                    
                                    r_list.append(input_mesh[i])
                                    psi_list.append(0.0)
                                }

                            }
                            else{

                                if(state_energy < 0.003){
                                    state_energy = 0.0
                                }
                                potential_schrod_eq_subroutine_inst.schroedinger_photoionization_subroutine(z_value: z_value, trial_energy: state_energy, number_blocks: number_of_blocks, l_number: i.quant_l-1, number_points: number_points, initial_delta_x: delta_x_initial, mesh_scalar: scalar, principal_quant_number: i.quant_n, input_pot_list: interpolated_Potential, thresh_criterion: thresh_criterion, max_thresh_iterations: max_thresh_iterations, left_energy_scalar: left_energy_scalar, right_energy_scalar: right_energy_scalar, r_list: interpolated_r_mesh)
                                
                                for i in potential_wavefunction_values_inst.norm_Pwavefunction_tuple_array {
                                    for j in i.wavefunction_list {
                                        r_list.append(j.r_value)
                                        psi_list.append(j.wavefunction_value)
                                    }
                                }

                            }
                            
                            }
                        //Unpack fine wavefunction to coarser grid for integration
                        if(r_list.count > 441){
//                            print("unpacking")
                            unpackedTuple = unpack_r_and_wavefunction(codex: interpolation_codex, r_values: r_list, psi_values: psi_list)
                            if(unpackedTuple.r_values.count == 0 || unpackedTuple.psi_values.count == 0){
                                unpackedTuple.r_values.append(0.0)
                                unpackedTuple.psi_values.append(0.0)
                            }
                        }
                        results_array_For_l_Minus_One.append((r_list: unpackedTuple.r_values, psi_list: unpackedTuple.psi_values, quant_n: i.quant_n, quant_l: i.quant_l, quant_m: i.quant_m, number_electrons: i.numb_electrons, photo_energy: photon_energy, new_energy: state_energy-photon_energy))
                        full_mesh = potential_mesh_potential_init_inst.r_mesh
                        potential_mesh_potential_init_inst.KE_term_list.removeAll()
                        potential_mesh_potential_init_inst.r_mesh = [0.0]
                        potential_mesh_potential_init_inst.trial_pot_list.removeAll()
                        potential_mesh_potential_init_inst.delta_r_list.removeAll()
                        potential_schrod_eq_subroutine_inst.test_final_values_list.removeAll()
                        potential_schrod_eq_subroutine_inst.test_final_r_values_list.removeAll()
                        potential_wavefunction_values_inst.outward_Pwavefunction_tuple_array.removeAll()
                        potential_wavefunction_values_inst.inward_Pwavefunction_tuple_array.removeAll()
                        potential_wavefunction_values_inst.inward_log_deriv_integral_values.removeAll()
                        potential_wavefunction_values_inst.outward_log_deriv_integral_values.removeAll()
                        potential_wavefunction_values_inst.norm_Pwavefunction_tuple_array.removeAll()

                        full_mesh = potential_mesh_potential_init_inst.r_mesh

                    }
                    
                    
                    
                }
                
//                // Calculates total electron density using wavefunction values and orbital values.  Using electron density, Z value, exchange alpha constant, potential values and radius values, modified potential is calculated
//                let rtot_array: [Double] = self.calculate_total_charge_density(total_mesh_count: total_mesh_count, core_shells: core_shells, val_shells: val_shells, results_array: results_array)
//                let pot_r_list: [Double] = full_mesh
//                var pot_final: [Double] = self.modified_potential(pot_r_list: pot_r_list, rtot_array: rtot_array, z_value: z_value, exchange_alpha: exchange_alpha, ionticity: ionticity)
//                
//                // If its an ion, uses a slightly different program to calculate final modified potential
//                if ionic_radius > 0.0 {
//                    pot_final = self.calculate_ionic_potential(pot_r_list: pot_r_list, pot_list: pot_final, ionticity: ionticity, branching_ratio: branching_ratio, ionic_radius: ionic_radius)
//                }
//                // Initializes and calculates beta values.  Finds maximum value from array
//                var beta_list: [Double] = [0.0]
//                for i in 0..<pot_r_list.count-1 {
//                    beta_list.append(abs(input_pot_list[i] - pot_final[i]))
//                }
//                current_beta_max = beta_list.max()!
//                
//                // Compares max calculated beta value to beta criterion.  If it's larger, modify initial potential using pratt scheme (in herman skillman book).
//                if (current_beta_max > beta_criterion) || (current_beta_max == 0.0) {
//                    
//                    for i in 0..<pot_r_list.count-1 {
//                        input_pot_list[i] = pratt_alpha*input_pot_list[i] + (1.0 - pratt_alpha)*pot_final[i]
//                    }
//                    
//                    print("CURRENT BETA: \(current_beta_max)")
//                    
//                    if current_beta_max == 0.0 {
//                        print(input_pot_list)
//                        for i in results_array {
//                            print(i.new_energy)
//                        }
//                    }
//                    
//                }
//                
//                beta_counter += 1
                
//            }
            
            print("")
            print("ITERATIONS: \(beta_counter) FINAL BETA: \(current_beta_max)")
            for i in results_array {
                print(i.new_energy)
            }
            print("")
            
        
//            print(results_array_For_l_Plus_One)
//            print(results_array_For_l_Minus_One)
            return(results_array_For_l_Plus_One, results_array_For_l_Minus_One, user_input_pot_list, input_mesh)
            
        }
    
}
