//
//  R Matrix Element Subroutine.swift
//  Atomic Structure Calculations (Non-Relativistic)
//
//  Created by Matthew Malaker on 4/19/22.
//

import Foundation
import Accelerate
class R_Element_Calculator: NSObject {
    
    
    ///calculateR
    ///Arguments:
    /// - Psi_I: The initial wavefunction in form of r, psi(r), n, l, m, numer_electorns, new_energy
    /// - Psi_F: The final wavefunction in form of r, psi(r), n, l, m, number_electrons, photon_energy, new_energy
    func calculateR(Psi_I:(r_list: [Double], psi_list: [Double], quant_n: Double, quant_l: Double, quant_m: Double, number_electrons: Double, new_energy: Double), Psi_F: (r_list: [Double], psi_list: [Double], quant_n: Double, quant_l: Double, quant_m: Double, number_electrons: Double, photo_energy: Double, new_energy: Double), plusOrMinus: Bool)->Double{

        var integral = 0.0
        var len = 1

        if(Psi_I.r_list.count < Psi_F.r_list.count){
            len = Psi_I.r_list.count-2
        }
        else{
            len = Psi_F.r_list.count-2
        }
        
        //Psis calculated are in hte form of r*psi(r), not (psi(r), so the r component needs to be divided out
        //This leads to NaNs at index 0, which should be zero in reality.
        //Accelerate is used to, well, accelerate this process
        
//        print(Psi_I.quant_n, Psi_I.quant_l, Psi_I.quant_m, Psi_F.photo_energy, Psi_F.psi_list.count, Psi_F.r_list.count)
        var initial_psi = vDSP.divide(Psi_I.psi_list, Psi_I.r_list)
        var final_psi = vDSP.divide(Psi_F.psi_list, Psi_F.r_list)
        initial_psi[0] = 0.0
        final_psi[0] = 0.0
        
        
        var fa = initial_psi[0]*Psi_I.r_list[0]*final_psi[0]
        for index in stride(from: 0, to: len, by: 1){
             
            let fb = initial_psi[index+1]*Psi_I.r_list[index+1]*final_psi[index+1]
            integral += 0.5*(Psi_I.r_list[index+1]-Psi_I.r_list[index])*(fa + fb)
            fa = fb // Assign right value fb to fa in order to avoid recaclulating it
        }
        
        return integral
    }
    
    func calculate_sigma(hv: Double, R_l_Plus_One: Double, R_l_Minus_One: Double, l: Double, num_Electrons: Double)->Double{
        
        //factor of two to convert hv from hartrees to rydbergs
        let alphanaught = 0.0072973525693
        let anaught = 5.29177210903*pow(10.0,-9.0)
        let const = (4.0*Double.pi*Double.pi*alphanaught*pow(anaught,2.0)/3.0)
        let value = const*hv*(num_Electrons/((2.0*l)+1.0))*(l*pow(R_l_Minus_One,2.0)+(l+1.0)*pow(R_l_Plus_One,2.0))*pow(10,18.0)
        return(value)
    }
    
    ///stacksigma
    ///
    ///Arguments:
    /// - insertStruct: An array of tuples containing the cross section for each energy for each nlm
    ///Returns:
    /// - stackedSigma: the sum of the cross section for each energy at each nlm. Is sigma(hv) totaled for all orbitals
    func stacksigma(insertStruct: [(n: Double, l: Double, m: Double, hv: [Double], sigma: [Double])])->(hv: [Double], sigma: [Double]){
        
        
        //We will not have a cross section at each hv for each nlm. This is because th photons can be too low of energy to ionize
        //Those cross sections are not even calculated, so there aren't even elements in insertStruct
        //We thus need to search insertstruct for values corresponding to each hv and add them together and
        //not get out-f-bounds errors
        
        let hvList = insertStruct.last!.hv
        
        var stackedSigma: (hv: [Double],sigma: [Double]) = (hv: [], sigma: [])
        var hvArray: [Double] = []
        var sigmaArray: [Double] = []
        
//        print(insertStruct)
        var sigma = 0.0
        for hv in hvList{
            for nlm in insertStruct{
                
                let index = nlm.hv.firstIndex(where: {$0 == hv})
                if(index ?? -1 >= 0){
                    sigma += nlm.sigma[index!]
                }
            }
            hvArray.append(hv)
            sigmaArray.append(sigma)
            sigma = 0.0
            
        }
        
        
//        for energy in stride(from: 0, to: insertStruct.last!.hv.count-1, by: 1){
////            let energy = insertStruct[0].hv[energy]
//            var sigma = 0.0
//            for index in stride(from: 0, to: insertStruct.count-1, by: 1){
//
//                sigma += insertStruct[index].sigma[energy]
//            }
//            hvArray.append(insertStruct[0].hv[energy])
//            sigmaArray.append(sigma)
//        }
        
        stackedSigma = (hv: hvArray, sigma: sigmaArray)
        return stackedSigma
    }
    
}
