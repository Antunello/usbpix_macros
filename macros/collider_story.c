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
#include "TLatex.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TLegend.h"
#include "THStack.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TColor.h"
#include <time.h>
#include "TGraphErrors.h"

void collider_story(std::string filename){
	ifstream dataFile;
	dataFile.open(filename.c_str());
	int counter = 0;
	std::string experiment;
	double start_dc, stop_dc, start_op, stop_op;
	std::vector<double> v_start_dc, v_stop_op;
	std::vector<double> v_counter, v_counter_err, v_mean_dc, v_err_dc, v_mean_op, v_err_op;
	std::vector<std::string> v_experiment;

	std::stringstream buffer;
	//std::cout<<dataFile.rdbuf();
	buffer << dataFile.rdbuf();
	while(buffer){
		counter++;
		std::cout<<buffer<<std::endl;
		buffer>>experiment;
		buffer>>start_dc;
		buffer>>stop_dc;
		buffer>>start_op;
		buffer>>stop_op;
			
		v_experiment.push_back(experiment);
		v_mean_dc.push_back((start_dc+stop_dc)/2);
		v_start_dc.push_back(start_dc);
		v_err_dc.push_back((stop_dc-start_dc)/2);
		v_mean_op.push_back((start_op+stop_op)/2);
		v_stop_op.push_back(stop_op);
		v_err_op.push_back((stop_op-start_op)/2);
		v_counter.push_back(counter);
		v_counter_err.push_back(0.3);
	}	
	std::reverse(v_counter.begin(), v_counter.end());
	TGraphErrors *op_graph = new TGraphErrors(v_counter.size()-1, &(v_mean_op[0]), &(v_counter[0]), &(v_err_op[0]), &(v_counter_err[0]));
	TGraphErrors *dc_graph = new TGraphErrors(v_counter.size()-1, &(v_mean_dc[0]), &(v_counter[0]), &(v_err_dc[0]), &(v_counter_err[0]));

	op_graph->SetFillColor(kPink-5);
	dc_graph->SetFillColor(kAzure+7);

	TMultiGraph *multi_graph = new TMultiGraph();
	multi_graph->Add(op_graph,"e2");
	multi_graph->Add(dc_graph, "e2");
	multi_graph->SetTitle(";year;");
	TCanvas *c1 = new TCanvas("c1","c1", 800, 600);
	multi_graph->Draw("AE2");
	multi_graph->GetYaxis()->SetLabelOffset(10);
	multi_graph->GetYaxis()->SetTickSize(0);
	multi_graph->Draw("AE2");

	TLegend *leg = new TLegend(0.7,0.72,0.9,0.87);	
	leg->AddEntry(dc_graph, "Construction" ,"F");
	leg->AddEntry(op_graph, "Operation" ,"F");
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->Draw();
  TLatex *lx1 = new TLatex();
  lx1->SetNDC(kFALSE);
  lx1->SetTextSize(0.04);
  lx1->SetTextAlign(12);
	lx1->SetTextColor(kBlack);
  //lx1->SetTextAngle(90);
	TLatex *lx2 =(TLatex*) lx1->Clone();
	lx2->SetTextColor(kWhite);
	
	for(int i=0; i<v_counter.size() -1; i++){

		//double time = (v_stop_op[i]+v_start_dc[i])/2 ;
		double time = v_start_dc[i] + 0.5 ;
		lx1->DrawText(time, v_counter[i],v_experiment[i].c_str());
		lx2->DrawText(time+0.1, v_counter[i]+0.01,v_experiment[i].c_str());
	}
	lx1->Draw();
	lx2->Draw();	
	gPad->SetGridx();
}

