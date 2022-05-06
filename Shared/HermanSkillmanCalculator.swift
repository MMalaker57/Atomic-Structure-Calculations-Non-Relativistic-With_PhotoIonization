//
//  HermanSkillmanCalculator.swift
//  HermanSkillmanCalculator
//
//  Created by Jeff_Terry on 11/29/21.
//

import Foundation
import SwiftUI
import CorePlot

class HermanSkillmanCalculator: ObservableObject {
    
    var plotDataModel: PlotDataClass? = nil
    // Class instances
    let functional_functions_inst = Functional_Functions()
    let hfs_textbook_values_inst = HFS_textbook_values()
    let wavefunction_values_inst = WaveFunction_Values()
    let schrod_eq_subroutine_inst = Schroedinger_Eq_SubRoutine()
    let mesh_potential_init_inst = Mesh_Potential_Initialization()
    let self_consistent_potentials_inst = Self_Consistent_Potentials()
    let r_calculator = R_Element_Calculator()
    let input_file_inst = Input_File()
    let output_file_inst = Output_File()
    var orbital_energies :[String] = []
    var element = ""
    var inputFileString = ""
    
    var alpha = 0.75
    
    let file_path_input = "file:///Users/jterry94/Downloads/input6.dat"
    
    // Arrays that store all calculated values including energy, potential and weave function values
    var results_array: [(r_list: [Double], psi_list: [Double], quant_n: Double, quant_l: Double, quant_m: Double, number_electrons: Double, new_energy: Double)] = []
    var results_array_For_l_Plus_One: [(r_list: [Double], psi_list: [Double], quant_n: Double, quant_l: Double, quant_m: Double, number_electrons: Double, photo_energy: Double, new_energy: Double)] = []
    var results_array_For_l_Minus_One: [(r_list: [Double], psi_list: [Double], quant_n: Double, quant_l: Double, quant_m: Double, number_electrons: Double, photo_energy: Double, new_energy: Double)] = []
    var results_array_For_l_Plus_One_Pos_Only: [(r_list: [Double], psi_list: [Double], quant_n: Double, quant_l: Double, quant_m: Double, number_electrons: Double, photo_energy: Double, new_energy: Double)] = []
    var results_array_For_l_Minus_One_Pos_Only: [(r_list: [Double], psi_list: [Double], quant_n: Double, quant_l: Double, quant_m: Double, number_electrons: Double, photo_energy: Double, new_energy: Double)] = []
    var sigma_results_struct: [(n: Double, l: Double, m: Double, hv: [Double], sigma: [Double])] = []
    var stackedSigma: (hv: [Double], sigma: [Double]) = (hv: [], sigma: [])
    var R_l_Plus_Or_Minus_One_Tuple_Array: [(n: Double, l: Double, m: Double, hv: Double, energy: Double, R_l_Plus_One: Double, R_l_Minus_One: Double, number_electrons: Double)] = []
    
    var sigma_array:[(n: Double, l: Double, Double, cross_sections: [Double])] = []
    
    var self_consistent_pot_list: [Double] = []
    var full_mesh: [Double] = []
    
    
    
    // Calculate NonRelativisticHartreeFock:
    // After parameters have been set and the input file has been selected, this function starts the Herman-Skillman program.  Messages and values will be printed in the terminal to help inform the user such as the beta value of the current iteration
    func calculateNonRelativisticHartreeFock() {
        
      
        
        //output: (atom name, Z, exchange_alpha, KEY, BETA, THRESH, mesh_count, MAXIT, KUT, RADION, BRATIO, XION, # core shells, #val shells, potential values, electron configuration)
        //let read_input_file_tuple = input_file_inst.read_DAT_110pt_inputfile(input_file_path: URL(string: file_path_input)!)
        
            
        let read_input_file_tuple = input_file_inst.read_DAT_110pt_inputfile(input_file_string: inputFileString)
        
        
        let atom_name: String = read_input_file_tuple.0
        let z_value: Double = read_input_file_tuple.1
        let exchange_alpha: Double = read_input_file_tuple.2
        let KEY: Int = read_input_file_tuple.3
        let beta_criterion: Double = read_input_file_tuple.4
        let thresh_criterion: Double = read_input_file_tuple.5
        let mesh_count: Double = read_input_file_tuple.6
        let max_beta_iterations: Int = read_input_file_tuple.7
        let ionic_radius: Double = read_input_file_tuple.9
        let branching_ratio: Double = read_input_file_tuple.10
        let ionticity: Double = read_input_file_tuple.11
        let core_shells: Int = read_input_file_tuple.12
        let val_shells: Int = read_input_file_tuple.13
        let user_input_pot_list: [Double] = read_input_file_tuple.14
        let electron_config_array: [(quant_n: Double, quant_l: Double, quant_m: Double, numb_electrons: Double, trial_energy: Double)] = read_input_file_tuple.15
        let photon_energy_array: [Double] = read_input_file_tuple.16
        
        let scalar: Double = (1.0/2.0)*pow((3.0*Double.pi)/4.0, 2.0/3.0)*pow(z_value, -1.0/3.0)
        let max_thresh_iterations: Int = 21
        let number_of_blocks: Double = (mesh_count - 1.0)/40.0
        let number_points: Double = 40.0
        let delta_x_initial: Double = 0.0025
        //var pratt_alpha: Double = 0.9
        var pratt_alpha: Double = 0.7
        
        
        
        element = atom_name
        
        pratt_alpha = alpha
        
        self_consistent_potentials_inst.myHartreeFockSCFCalculator = self
        output_file_inst.myHartreeFockSCFCalculator = self

        
        // Times the main section of the code starting from this line
        let t1 = mach_absolute_time()
        
        // Array tuple that calls main program and stores its values
        let results_array_tuple = self_consistent_potentials_inst.self_consistent_potential(z_value: z_value, delta_x_initial: delta_x_initial, number_of_blocks: number_of_blocks, scalar: scalar, number_points: number_points, user_input_pot_list: user_input_pot_list, core_shells: core_shells, val_shells: val_shells, electron_config_array: electron_config_array, beta_criterion: beta_criterion, exchange_alpha: exchange_alpha, pratt_alpha: pratt_alpha, KEY: KEY, thresh_criterion: thresh_criterion, max_beta_iterations: max_beta_iterations, max_thresh_iterations: max_thresh_iterations, ionic_radius: ionic_radius, branching_ratio: branching_ratio, ionticity: ionticity, photoionization_energies: photon_energy_array)
        
        print(electron_config_array)
        var electron_config_array_2 = electron_config_array
        for i in stride(from: 0, to: results_array_tuple.0.count, by: 1){
//            print(results_array_tuple.0[i].new_energy)
            electron_config_array_2[i].trial_energy = results_array_tuple.0[i].new_energy

        }
        
        results_array = results_array_tuple.0
        self_consistent_pot_list = results_array_tuple.1
        full_mesh = results_array_tuple.2
        
        let results_array_tuple_photo = self_consistent_potentials_inst.self_consistent_potential2(z_value: z_value, delta_x_initial: delta_x_initial, number_of_blocks: number_of_blocks, scalar: scalar, number_points: number_points, user_input_pot_list: self_consistent_pot_list, core_shells: core_shells, val_shells: val_shells, electron_config_array: electron_config_array_2, beta_criterion: beta_criterion, exchange_alpha: exchange_alpha, pratt_alpha: pratt_alpha, KEY: KEY, thresh_criterion: thresh_criterion, max_beta_iterations: max_beta_iterations, max_thresh_iterations: max_thresh_iterations, ionic_radius: ionic_radius, branching_ratio: branching_ratio, ionticity: ionticity, photoionization_energies: photon_energy_array, r_mesh: full_mesh)
        
        // Separates 'results_array_tuple' into five arrays, wavefunction/energy values, potential values, r mesh values, and l+-1 photoionization wavefunctions
        results_array = results_array_tuple.0
        self_consistent_pot_list = results_array_tuple.1
        full_mesh = results_array_tuple.2
        
        results_array_For_l_Plus_One = results_array_tuple_photo.0
        results_array_For_l_Minus_One = results_array_tuple_photo.1
        
        print(results_array_For_l_Minus_One.count, results_array_For_l_Plus_One.count)
        
        
        // Ends timing of main section of code
        let t2 = mach_absolute_time()
        
        // Calculates and prints run time of the main code
        let elapsed = t2 - t1
        var timeBaseInfo = mach_timebase_info_data_t()
        mach_timebase_info(&timeBaseInfo)
        let elapsedNano = elapsed * UInt64(timeBaseInfo.numer) / UInt64(timeBaseInfo.denom);
        
        print("Execution time for the code is:  \(Double(elapsedNano)*1.0E-09) seconds")
        print("")
        
        
        
        //Calculates R+-1 for each energy in each nlm to feed into photoionization cross section
        for i in stride(from: 0, through: results_array_For_l_Plus_One.count-1, by: 1){
            
            let n = results_array_For_l_Plus_One[i].quant_n
            let l = results_array_For_l_Plus_One[i].quant_l
            let m = results_array_For_l_Plus_One[i].quant_m
            let photon_energy = results_array_For_l_Plus_One[i].photo_energy
            
            //check for which initial psi has the same nlm as the final
            //nlm is used as a signature. We do not incorporate the l+-1 in that tag
            var results_Array_Index = 0 // Needs to be Int
            for element in stride(from:0, through:results_array.count-1, by: 1){
                if results_array[element].quant_n == n && results_array[element].quant_l == l && results_array[element].quant_m == m{
                    results_Array_Index = element
                }
            }
            
            
            if(results_array_For_l_Minus_One[i].r_list.count == 0){
                results_array_For_l_Minus_One[i].r_list.append(0.0)
                results_array_For_l_Minus_One[i].psi_list.append(0.0)
            }
//            print("Initial Array is \(results_array.count), final + 1 is \(results_array_For_l_Plus_One.count) long, and final - 1 is \(results_array_For_l_Minus_One.count) long.")
            
            //Might have to change the energy portion of this tuple
            R_l_Plus_Or_Minus_One_Tuple_Array.append((n: n, l: l, m: m, hv: results_array_For_l_Plus_One[i].photo_energy, energy: results_array_For_l_Plus_One[i].new_energy, R_l_Plus_One: r_calculator.calculateR(Psi_I: results_array[results_Array_Index], Psi_F: results_array_For_l_Plus_One[i], plusOrMinus: true), R_l_Minus_One: r_calculator.calculateR(Psi_I: results_array[results_Array_Index], Psi_F: results_array_For_l_Minus_One[i], plusOrMinus: false), number_electrons: results_array_For_l_Minus_One[i].number_electrons))
            
            
        }
        
        R_l_Plus_Or_Minus_One_Tuple_Array = R_l_Plus_Or_Minus_One_Tuple_Array.sorted(by: {($0.n, $0.l, $0.m, $0.hv) < ($1.n, $1.l, $1.m, $1.hv)})
//        print(R_l_Plus_Or_Minus_One_Tuple_Array)
        
        //calculates sigma for each energy at each nlm

        var sigma_n = R_l_Plus_Or_Minus_One_Tuple_Array[0].n
        var sigma_l = R_l_Plus_Or_Minus_One_Tuple_Array[0].l
        var sigma_m = R_l_Plus_Or_Minus_One_Tuple_Array[0].m
        
        var num_elec = R_l_Plus_Or_Minus_One_Tuple_Array[0].number_electrons
        var sigma_array: [Double] = []
        var hv_array: [Double] = []
        
        for sigma in stride(from: 0, to:R_l_Plus_Or_Minus_One_Tuple_Array.count-1,by: 1){

//            print("sigma_n = \(sigma_n)")
            var hv = R_l_Plus_Or_Minus_One_Tuple_Array[sigma].hv

            let tempSigma = (r_calculator.calculate_sigma(hv: hv, R_l_Plus_One: R_l_Plus_Or_Minus_One_Tuple_Array[sigma].R_l_Plus_One, R_l_Minus_One: R_l_Plus_Or_Minus_One_Tuple_Array[sigma].R_l_Minus_One, l: sigma_l, num_Electrons: num_elec))
            
            if(tempSigma>0.0 /*&& tempSigma.isNaN == false*/){
            sigma_array.append(tempSigma)
                hv_array.append(hv)
            }

            if(R_l_Plus_Or_Minus_One_Tuple_Array[sigma+1].n != sigma_n || R_l_Plus_Or_Minus_One_Tuple_Array[sigma+1].l != sigma_l || R_l_Plus_Or_Minus_One_Tuple_Array[sigma+1].m != sigma_m || sigma == R_l_Plus_Or_Minus_One_Tuple_Array.count-2){
                
                if(hv_array != []){
                    sigma_results_struct.append((n: sigma_n, l: sigma_l, m: sigma_m, hv: hv_array, sigma: sigma_array))
                }
                sigma_n = R_l_Plus_Or_Minus_One_Tuple_Array[sigma+1].n
                sigma_l = R_l_Plus_Or_Minus_One_Tuple_Array[sigma+1].l
                sigma_m = R_l_Plus_Or_Minus_One_Tuple_Array[sigma+1].m
                sigma_array = []
                hv_array = []
                num_elec = R_l_Plus_Or_Minus_One_Tuple_Array[sigma+1].number_electrons
            }
        }

        // Prints r mesh and corresponding potential value
        for i in 0..<self_consistent_pot_list.count {
            print(full_mesh[i], full_mesh[i]/scalar, self_consistent_pot_list[i]/(-2.0*z_value))
        }
        
        // Removes all values from 'orbital_energies' drop down menu
        orbital_energies.removeAll()
        
        for i in results_array_For_l_Plus_One{
//            print("n=\(i.quant_n), l = \(i.quant_l), state_energy = \(i.photo_energy + i.new_energy)")
            if(i.photo_energy+i.new_energy) >= 0.0{
//                print("triggered addition to positive only")
                results_array_For_l_Plus_One_Pos_Only.append(i)
            }
        }
        
        for i in results_array_For_l_Minus_One{
            if(i.photo_energy+i.new_energy) >= 0.0{
                results_array_For_l_Minus_One_Pos_Only.append(i)
            }
        }
        
        
        
        // Adds newly calculated orbital energies with correesponding quantum numbers as items in drop down menu.
        for i in results_array {
            orbital_energies.append("\(Int(i.quant_n))\(Int(i.quant_l))\(Int(i.quant_m)): \(i.new_energy)")
        }
//        print(sigma_results_struct)
        for i in results_array_For_l_Plus_One_Pos_Only{
            var kinetic_energy = i.photo_energy + i.new_energy
            if(kinetic_energy < 0.0012){
                kinetic_energy = 0.0
            }
                orbital_energies.append("\(Int(i.quant_n))\(Int(i.quant_l+1.0))\(Int(i.quant_m)): Final State Wf Ekin = \(kinetic_energy)")
            
        }
        
        for i in sigma_results_struct{
            orbital_energies.append("\(Int(i.n))\(Int(i.l))\(Int(i.m)): PhotoIonization")
            
        }
        orbital_energies.append("Stacked PhotoIonization")
        
        stackedSigma = r_calculator.stacksigma(insertStruct: sigma_results_struct)
        
        // Creates an energy output file and wave function/potential output file.
        output_file_inst.make_energy_output_file(atom_name: atom_name, z_value: z_value, final_results_array: results_array)
        output_file_inst.photoIonOutput(atom_name: atom_name, z_value: z_value, photoIon_array: sigma_results_struct, stacked_sigma: stackedSigma)
        output_file_inst.make_wf_pot_output_file(atom_name: atom_name, z_value: z_value, final_results_array: results_array, self_consistent_pot_list: self_consistent_pot_list, r_mesh: full_mesh)
        
    }
}
