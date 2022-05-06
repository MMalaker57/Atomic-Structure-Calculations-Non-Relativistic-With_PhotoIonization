//
//  ContentView.swift
//  Shared
//
//  Created by Jeff Terry on 11/28/21.
//

import SwiftUI
import CorePlot

typealias plotDataType = [CPTScatterPlotField : Double]

struct ContentView: View {
    
    @ObservedObject private var dataCalculator = CalculatePlotData()
    @EnvironmentObject var plotDataModel :PlotDataClass
    @ObservedObject var hartreeFockSCFCalculator = HermanSkillmanCalculator()
    
    @State var selectedOutputIndex :Int = 0
    @State var selectedWavefunctionIndex :Int = -1
    
    @State var selectedWavefunction = ""
    @State var filename = ""
    
    @State var outputPickerArray :[String] = []
    @State var outputWavefunctionArray :[String] = []
    
    @State var isImporting: Bool = false
    
    @State var element = ""
    
    @State var alphaParameter = "0.75"
    @State var inputFileString = ""
    
    @State var potential :[Double] = []
    @State var wavefunctionResults: [(r_list: [Double], psi_list: [Double], quant_n: Double, quant_l: Double, quant_m: Double, number_electrons: Double, new_energy: Double)] = []
    @State var wavefunctionResults_lPlusOne: [(r_list: [Double], psi_list: [Double], quant_n: Double, quant_l: Double, quant_m: Double, number_electrons: Double, photo_energy: Double, new_energy: Double)] = []
    @State var wavefunctionResults_lMinusOne: [(r_list: [Double], psi_list: [Double], quant_n: Double, quant_l: Double, quant_m: Double, number_electrons: Double, photo_energy: Double, new_energy: Double)] = []
    @State var sigma_results: [(n: Double, l: Double, m: Double, hv: [Double], sigma: [Double])] = []
    @State var stacked_results: (hv: [Double], sigma: [Double]) = (hv: [], sigma: [])
    @State var mesh: [Double] = []
    

    var body: some View {
        
        VStack{
      
            CorePlot(dataForPlot: $plotDataModel.plotData, changingPlotParameters: $plotDataModel.changingPlotParameters)
                .setPlotPadding(left: 10)
                .setPlotPadding(right: 10)
                .setPlotPadding(top: 10)
                .setPlotPadding(bottom: 10)
                .padding()
            
            
            
            Divider()
            
            HStack{
                
                Button("Load Input File", action: {
                                isImporting = false
                                
                                //fix broken picker sheet
                                DispatchQueue.main.asyncAfter(deadline: .now() + 0.1) {
                                    isImporting = true
                                }
                            })
                                .padding()
                                .fileImporter(
                                            isPresented: $isImporting,
                                            allowedContentTypes: [.text],
                                            allowsMultipleSelection: false
                                        ) { result in
                                            do {
                                                guard let selectedFile: URL = try result.get().first else { return }
                                                
                                                print("Selected file is", selectedFile)
                                                
                                                //trying to get access to url contents
                                                if (CFURLStartAccessingSecurityScopedResource(selectedFile as CFURL)) {
                                                                        
                                                    filename = selectedFile.lastPathComponent
                                                    
                                                    guard let message = String(data: try Data(contentsOf: selectedFile), encoding: .utf8) else { return }
                                                        
                                                    //print(message)
                                                    
                                                    inputFileString = message
                                                        
                                                    //done accessing the url
                                                    CFURLStopAccessingSecurityScopedResource(selectedFile as CFURL)
                                                    
                                                    
                                                }
                                                else {
                                                    print("Permission error!")
                                                }
                                            } catch {
                                                // Handle failure.
                                                print(error.localizedDescription)
                                            }
                                        }
                
                Text($filename.wrappedValue)
                
                
            }

            
            Divider()

            
            
            HStack{
                
                VStack{
                    
                    Text("Alpha Parameter")
                    
                    TextField("Alpha Parameter between 0.1 and 0.95", text: $alphaParameter)
                        .frame(minWidth: 100.0, idealWidth: 300.0, maxWidth: 400.0, alignment: .center)
                }
                
                Button("Calculate", action: {self.calculateHartreeFock()} )
                .padding()
            
            }
            
            Divider()
            
            HStack{
                
               // Spacer()
                
#if os(macOS)
                
                Picker("Output", selection: $selectedWavefunction, content: {
                    ForEach(outputWavefunctionArray, id:\.self) {
                    //ForEach(0..<outputWavefunctionArray.count) {
                        Text($0).tag($0)
                                        }
                })
                    .padding()
                    .frame(minWidth: 100.0, idealWidth: 300.0, maxWidth: 400.0, alignment: .center)
                    .onChange(of: selectedWavefunction, perform: { selection in
                        
                        if let index = outputWavefunctionArray.firstIndex(of: selection) {
                            
                            selectedWavefunctionIndex = index
                            
                        }
                        
                        if selectedWavefunctionIndex < (outputWavefunctionArray.count - 2){
                            
                            updateWavefunctionPlot(index: selectedWavefunctionIndex)
                            
                        }
                        else{
                            if selectedWavefunctionIndex == (outputWavefunctionArray.count - 2){
                                updateStackedPlot()
                            }
                            else{
                            updatePotentialPlot()
                            }
                            
                        }
                    })
    
#elseif os(iOS)
                
                NavigationView{
                    
                    Form{
                        
                        Section{
                            
                            Picker("Output", selection: $selectedWavefunction, content: {
                                ForEach(outputWavefunctionArray, id:\.self) {
                                //ForEach(0..<outputWavefunctionArray.count) {
                                    Text($0).tag($0)
                                                    }
                            })
                                .padding()
                                .frame(minWidth: 100.0, idealWidth: 300.0, maxWidth: 400.0, alignment: .center)
                                .onChange(of: selectedWavefunction, perform: { selection in
                                    
                                    if let index = outputWavefunctionArray.firstIndex(of: selection) {
                                        
                                        selectedWavefunctionIndex = index
                                        print("index = \(index)")
                                        
                                    }
                                    
                                    if selectedWavefunctionIndex < (outputWavefunctionArray.count - 2){
                                        
                                        updateWavefunctionPlot(index: selectedWavefunctionIndex)
                                        print("index = \(selectedWavefunctionIndex) and ")
                                        
                                    }
                                    else{
                                        print("index = \(selectedWavefunctionIndex)")
                                        if(selectedWavefunctionIndex == outputWavefunctionArray.count - 2){
                                            updateStackedPlot()
                                            
                                        }
                                        if(selectedWavefunctionIndex == outputWavefunctionArray.count-1){
                                            updatePotentialPlot()
                                        }
                                    }
                                })
                        }
                    }
                }
#endif
            }
            
            Spacer()
        }
    }
    
    func updateStackedPlot(){
        print(stacked_results)
        if mesh.count == 0 {
            
            return
        }
        
        //set the Plot Parameters
        
        plotDataModel.changingPlotParameters.xLabel = "hv"
        plotDataModel.changingPlotParameters.yLabel = "sigma"
        plotDataModel.changingPlotParameters.lineColor = .red()
        plotDataModel.changingPlotParameters.title = "Stacked Photoionization"
        
        plotDataModel.zeroData()
    
        let MaxXDisplayCoord = stacked_results.hv.max()!
        let MinXDisplayCoord = stacked_results.hv.min()!
        let MinYDisplayCoord = stacked_results.sigma.min()! - 1
        
        let MaxYDisplayCoord = stacked_results.sigma.max()! * 1.25
        
        plotDataModel.changingPlotParameters.xMax = MaxXDisplayCoord
        plotDataModel.changingPlotParameters.yMax = MaxYDisplayCoord
        plotDataModel.changingPlotParameters.yMin = MinYDisplayCoord
        plotDataModel.changingPlotParameters.xMin = MinXDisplayCoord
//        print("sigma_results.count = \(sigma_results.count)")
        
        for i in stride(from: 0, through: stacked_results.hv.count-1, by: 1) {
            let x = stacked_results.hv[i]
            let y = stacked_results.sigma[i]
            if abs(y) < 50.0{
                let dataPoint: plotDataType = [.X: x, .Y: y]
                plotDataModel.plotData.append(dataPoint)
            }
        }
        
    }
    
    func updatePotentialPlot(){
        print(selectedWavefunctionIndex)
        if mesh.count == 0 {
            
            return
        }
        
    
        
        //set the Plot Parameters
        plotDataModel.changingPlotParameters.xMin = -3.0
        plotDataModel.changingPlotParameters.xLabel = "r"
        plotDataModel.changingPlotParameters.yLabel = "V(r)"
        plotDataModel.changingPlotParameters.lineColor = .blue()
        plotDataModel.changingPlotParameters.title = element + " Potential"
        
        plotDataModel.zeroData()
        
        let MaxXDisplayCoord = mesh.max()!*0.75
        let MinYDisplayCoord = potential.min()!*1.1
        var MaxYDisplayCoord = potential.max()!*1.1
        
        if MaxYDisplayCoord < -0.0 {
            
            MaxYDisplayCoord = -1.0*MinYDisplayCoord/10.0
        }
        
        plotDataModel.changingPlotParameters.xMax = MaxXDisplayCoord
        plotDataModel.changingPlotParameters.yMax = MaxYDisplayCoord
        plotDataModel.changingPlotParameters.yMin = MinYDisplayCoord
        
        for i in stride(from: 0, to: mesh.count, by: 1) {
            
            let x = mesh[i]
            let y = potential[i]
            
            let dataPoint: plotDataType = [.X: x, .Y: y]
            plotDataModel.plotData.append(dataPoint)
        }
        
        
        
        
    }
    
    
    func updateWavefunctionPlot(index: Int){
        print(index, (wavefunctionResults.count))
        if wavefunctionResults.count == 0 {
            return
        }
        if(index < (wavefunctionResults.count)){
        print("triggered wf plotting")
        //set the Plot Parameters
        
        plotDataModel.changingPlotParameters.xLabel = "r"
        plotDataModel.changingPlotParameters.yLabel = "r 𝜳(r)"
        plotDataModel.changingPlotParameters.lineColor = .red()
        plotDataModel.changingPlotParameters.title = element + " " + outputWavefunctionArray[index] + " Ry"
        
        plotDataModel.zeroData()
    
        
        let MaxXDisplayCoord = wavefunctionResults[index].r_list.max()!*1.25
        var MinYDisplayCoord = wavefunctionResults[index].psi_list.min()!*1.25
        if MinYDisplayCoord > -0.11 {
            
            MinYDisplayCoord = -0.35
        }
        
        let MaxYDisplayCoord = wavefunctionResults[index].psi_list.max()!*1.25
        
        plotDataModel.changingPlotParameters.xMax = MaxXDisplayCoord
        plotDataModel.changingPlotParameters.yMax = MaxYDisplayCoord
        plotDataModel.changingPlotParameters.yMin = MinYDisplayCoord
        plotDataModel.changingPlotParameters.xMin = -1.0*MaxXDisplayCoord/10.0
        
        
        for i in stride(from: 0, to: wavefunctionResults[index].r_list.count, by: 1) {
            
            let x = wavefunctionResults[index].r_list[i]
            let y = wavefunctionResults[index].psi_list[i]
            
            if abs(y) < 50.0{
                let dataPoint: plotDataType = [.X: x, .Y: y]
                plotDataModel.plotData.append(dataPoint)
            }
        }
    
    }
        
        else{
            
            if index < wavefunctionResults.count + wavefunctionResults_lPlusOne.count{
                print("triggered l+1 wf plotting")
                print(wavefunctionResults_lPlusOne[index - wavefunctionResults.count].r_list.max(), wavefunctionResults_lPlusOne[index - wavefunctionResults.count].psi_list)
                plotDataModel.changingPlotParameters.xLabel = "r"
                plotDataModel.changingPlotParameters.yLabel = "r 𝜳(r)"
                plotDataModel.changingPlotParameters.lineColor = .red()
                plotDataModel.changingPlotParameters.title = element + " " + outputWavefunctionArray[index] + " Ry"
                
                plotDataModel.zeroData()
//                print(wavefunctionResults_lPlusOne[index - wavefunctionResults.count])
                
                let MaxXDisplayCoord = wavefunctionResults_lPlusOne[index - wavefunctionResults.count].r_list.max()!*1.25
                var MinYDisplayCoord = wavefunctionResults_lPlusOne[index - wavefunctionResults.count].psi_list.min()!*1.25
                if MinYDisplayCoord > -0.11 {
                    
                    MinYDisplayCoord = -0.35
                }
                
                let MaxYDisplayCoord = wavefunctionResults_lPlusOne[index - wavefunctionResults.count].psi_list.max()!*1.25
         
                
                
                //CHANGE THIS BACK
                plotDataModel.changingPlotParameters.xMax = 18
                plotDataModel.changingPlotParameters.yMax = MaxYDisplayCoord
                plotDataModel.changingPlotParameters.yMin = MinYDisplayCoord
                plotDataModel.changingPlotParameters.xMin = -1.0*MaxXDisplayCoord/10.0
                
                
                for i in stride(from: 0, to: wavefunctionResults_lPlusOne[index - wavefunctionResults.count].r_list.count, by: 1) {
                    
                    let x = wavefunctionResults_lPlusOne[index - wavefunctionResults.count].r_list[i]
                    let y = wavefunctionResults_lPlusOne[index - wavefunctionResults.count].psi_list[i]
                    
                    if abs(y) < 50.0{
                        let dataPoint: plotDataType = [.X: x, .Y: y]
                        plotDataModel.plotData.append(dataPoint)
                    }
                }
            
            }
            
            
            else{
                print("triggered sigma plotting")
//            print(index - wavefunctionResults.count)
            //set the Plot Parameters
            
            plotDataModel.changingPlotParameters.xLabel = "hv"
            plotDataModel.changingPlotParameters.yLabel = "sigma"
            plotDataModel.changingPlotParameters.lineColor = .red()
            plotDataModel.changingPlotParameters.title = element + " " + outputWavefunctionArray[index - wavefunctionResults.count - wavefunctionResults_lPlusOne.count] + " Photoionization"
            
            plotDataModel.zeroData()
        
//            print(sigma_results[index - wavefunctionResults.count])
                let MaxXDisplayCoord = sigma_results[index - wavefunctionResults.count - wavefunctionResults_lPlusOne.count].hv.max()!
                let MinXDisplayCoord = sigma_results[index - wavefunctionResults.count - wavefunctionResults_lPlusOne.count].hv.min()!
            let MinYDisplayCoord = sigma_results[index - wavefunctionResults.count - wavefunctionResults_lPlusOne.count].sigma.min()! * 1.25
            
            let MaxYDisplayCoord = sigma_results[index - wavefunctionResults.count - wavefunctionResults_lPlusOne.count].sigma.max()! * 1.25
            
            plotDataModel.changingPlotParameters.xMax = MaxXDisplayCoord
            plotDataModel.changingPlotParameters.yMax = MaxYDisplayCoord
            plotDataModel.changingPlotParameters.yMin = MinYDisplayCoord
            plotDataModel.changingPlotParameters.xMin = MinXDisplayCoord
//            print("sigma_results.count = \(sigma_results.count)")
            print(sigma_results[index - wavefunctionResults.count - wavefunctionResults_lPlusOne.count])
            
            for i in stride(from: 0, through: sigma_results[index - wavefunctionResults.count - wavefunctionResults_lPlusOne.count].hv.count-1, by: 1) {
//                print(i)
                
                let x = sigma_results[index - wavefunctionResults.count - wavefunctionResults_lPlusOne.count].hv[i]
                let y = sigma_results[index - wavefunctionResults.count - wavefunctionResults_lPlusOne.count].sigma[i]
                if abs(y) < 50.0{
                    let dataPoint: plotDataType = [.X: x, .Y: y]
                    plotDataModel.plotData.append(dataPoint)
                }
            }
        
        }
        }
    }
    
    func calculateHartreeFock(){
        
        outputPickerArray = []
        outputWavefunctionArray = []
        
        wavefunctionResults.removeAll()
        potential.removeAll()
        mesh.removeAll()
        
        plotDataModel.calculatedText = ""
        
        
        if inputFileString != "" {
            
            hartreeFockSCFCalculator.inputFileString = inputFileString
            
        }
        else{
            
            return
        }
        
        
        hartreeFockSCFCalculator.alpha = Double(alphaParameter)!
        
        //pass the plotDataModel to the calculator
        hartreeFockSCFCalculator.plotDataModel = self.plotDataModel
        //Calculate the new plotting data and place in the plotDataModel
        hartreeFockSCFCalculator.calculateNonRelativisticHartreeFock()
        
        outputWavefunctionArray = hartreeFockSCFCalculator.orbital_energies
        outputWavefunctionArray.append("Potential")
        print("outputarray is \(outputWavefunctionArray.count) long")
        
        print("trigger index is \(outputWavefunctionArray.count-3)")
        
        wavefunctionResults = hartreeFockSCFCalculator.results_array
        
        wavefunctionResults_lPlusOne = hartreeFockSCFCalculator.results_array_For_l_Plus_One_Pos_Only
//        print(wavefunctionResults_lPlusOne)
        
        potential = hartreeFockSCFCalculator.self_consistent_pot_list
        sigma_results = hartreeFockSCFCalculator.sigma_results_struct
        stacked_results = hartreeFockSCFCalculator.stackedSigma
        
        mesh = hartreeFockSCFCalculator.full_mesh
        
        element = hartreeFockSCFCalculator.element
        plotDataModel.fileName = element
//        
//        for i in stride(from: 0, through: mesh.count - 2, by: 1){
//            print(i, mesh[i+1]-mesh[i], potential[i+1]-potential[i])
//            
//        }
//
        
        
        
        updateWavefunctionPlot(index: 0)
    
    }
    
}

struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView()
    }
}

