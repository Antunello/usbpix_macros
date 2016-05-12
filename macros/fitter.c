#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "TGraph.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include <string>
#include <sstream>
#include <fstream>
#include "TMultiGraph.h"
#include "TPad.h"
#include "TLegend.h"
#include "THStack.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TColor.h"
#include <time.h>

void fitter(std::string data_filename){

	std::vector<double> energy_vector, dEdx_vector;
	double energy, dEdx, trash;
	ifstream dataFile;
	dataFile.open(data_filename.c_str());
	std::stringstream buffer;
	buffer << dataFile.rdbuf();
	while(buffer){
		buffer>>energy;
		buffer>>trash;
		buffer>>trash;
		buffer>>dEdx;
		//buffer>>trash;
		//buffer>>trash;
		//buffer>>trash;

		energy_vector.push_back(energy);
		dEdx_vector.push_back(dEdx);
	} 
	dataFile.close();
	TGraph *conversion = new TGraph(energy_vector.size(), &(energy_vector[0]), &(dEdx_vector[0]));
	conversion->SetTitle(";Energy;dE/dx");
	TCanvas *c1 = new TCanvas("c1","c1",800,600);
	conversion->SetMarkerColor(kBlack);
	conversion->SetMarkerSize(2);
	conversion->SetMarkerStyle(24);
	conversion->Draw("AP");
	conversion->Fit("pol4");		


}
