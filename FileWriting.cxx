/////////////////////////////////////////////////////////////////////
////															 ////
//// Example 1: Get a threshold map and plot it in 2D histogram  ////
//// Meant to provide examples on very low level how to get data ////
//// from RootDBD data files as used in USBpix.					 ////
//// 															 ////
//// Author: Marcp Rimoldi	(with hug support from Antonello ;-), thanks! )  ////
////															 ////
/////////////////////////////////////////////////////////////////////


#include <DataContainer/PixDBData.h>
#include "PixConfDBInterface/PixConfDBInterface.h"
#include "PixConfDBInterface/RootDB.h"
#include "PixController/PixScan.h"
#include "GeneralDBfunctions.h"

#include <TH2F.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include <TApplication.h>
#include <TStyle.h>

#include <math.h>
#include <limits>
#include <string>

#include <iomanip>

#include <iostream>
#include <fstream>
#include <sstream>

//root headers
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TEfficiency.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TH2F.h"
#include "THStack.h"
#include "TF1.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include <time.h>
//#include <mytime.h>

using namespace PixLib;
using namespace std;

string ScanName, GroupName, ModuleName;
bool scan_found = false;
double TIDValue(double CT){

	return CT* 1e5 * 4.7474;
}
void listScans(string fname) {	// this function explains how to loop over DBInquires and get the scan labels. Searching for correct one and plotting it.
	string sname, gname, mname;
	RootDB *db = new RootDB(fname.c_str());
	DBInquire *root = db->readRootRecord(1);
	for (recordIterator i = root->recordBegin(); i != root->recordEnd(); i++) {
		if ((*i)->getName() == "PixScanResult") {
			sname = (*i)->getDecName();
			getDecNameCore(sname);
			//cout << "NEW SCAN sname: \"" << sname << "\"" << endl;
			for (recordIterator ii = (*i)->recordBegin(); ii != (*i)->recordEnd(); ii++) {
				if ((*ii)->getName() == "PixModuleGroup") {
					gname = (*ii)->getDecName();
					getDecNameCore(gname);
					//cout << "  with module group gname: " <<  gname << endl;
					for (recordIterator iii = (*ii)->recordBegin(); iii != (*ii)->recordEnd(); iii++) {
						if ((*iii)->getName() == "PixModule") {
							mname = (*iii)->getDecName();
							getDecNameCore(mname);
				//			cout << "     with module mname: " <<  mname << endl;
						}
					}
				}
			}
		}

		if (sname == "DIGITAL_TEST"/*"fdac tune FdacVbn=30"*/)
		{
			ScanName = sname;
			GroupName = gname;
			ModuleName = mname;
			//cout << "found scan: " << ScanName << ":" << GroupName << ":" << ModuleName << endl;
			break;
		}
	}
	delete db;
}


void fillData(string fname, string &my_time, vector<double> &v, vector<string> &v_name) {
	string sname, gname, mname;
	scan_found = false;
	RootDB *db = new RootDB(fname.c_str());
	DBInquire *root = db->readRootRecord(1);
	for (recordIterator i = root->recordBegin(); i != root->recordEnd(); i++) {
		if ((*i)->getName() == "PixScanResult") {
			sname = (*i)->getDecName();
			getDecNameCore(sname);
			//cout << "NEW SCAN: \"" << sname << "\"" << endl;
			for (recordIterator ii = (*i)->recordBegin(); ii != (*i)->recordEnd(); ii++) {
				if ((*ii)->getName() == "PixModuleGroup") {
					gname = (*ii)->getDecName();
					getDecNameCore(gname);
					//cout << "  with module group: " <<  gname << endl;
					for (recordIterator iii = (*ii)->recordBegin(); iii != (*ii)->recordEnd(); iii++) {
						if ((*iii)->getName() == "PixModule") {
							mname = (*iii)->getDecName();
							getDecNameCore(mname);
							//cout << "     with module: " <<  mname << endl;
						}
					}
				}
			}
		}

		if ((*i)->getName() == "DCS_readings") {
			sname = (*i)->getDecName();
			getDecNameCore(sname);
			//cout << sname << endl;
			if (sname == "Timestamps")
			{
				//cout << sname << endl;
				for (fieldIterator j = (*i)->fieldBegin(); j != (*i)->fieldEnd(); j++)
				{
					//cout << (*j)->getName() << endl;
					string bla = (*j)->getName();
					//bla = bla.substr(9);
					//cout << bla << endl;
					if (bla == "get time")
					{
						//cout << rep << endl;
						string value = "";
						db->DBProcess((*j), READ, value);
						//cout << "got time stamp: " << value << endl;
						my_time = value;
					}
				}
			}
			if (sname == "DCS_readings")
			{
				int a = 0;
				//cout << sname << endl;
				for (fieldIterator j = (*i)->fieldBegin(); j != (*i)->fieldEnd(); j++)
				{
					string blubb = (*j)->getName();
					float value = 0;
					db->DBProcess((*j), READ, value);
					//cout << "DCS READING: " << blubb << " " << value << endl;
					v.push_back(value);
					v_name.push_back(blubb);
				}
			}
		}
	}
}
/*
void plotFDACd(string fname,vector<double> &good) { // how to get the data and plot it. Still to do: make histogram beatiful ;-)
	PixDBData data("scan", " ", " ");
	TH2F * myData = data.GetMap(0, PixScan::OCCUPANCY, -1);
	cout << endl << "got data" << endl << endl;
	TH2F * myPlot = new TH2F(*(myData);)//((TH2F*)data.GetMap(0, PixScan::SCURVE_MEAN,0))
	cout << endl << "created TH2F" << endl << endl;

	int numgood=0;
	for (int numpix=0; numpix < myPlot->GetSize();numpix++){
		if((myPlot->GetBinContent(numpix))==200) numgood++;
	}
	good.push_back(numgood);

	TApplication app("app", NULL, NULL);
	myPlot->Draw("COLZ");
	cout << "Press ctrl-C to stop" << endl;
	app.Run();
}
*/
void plotANALOGDIGITAL(string fname, string sname, string gname, string mname, vector<double> &good, vector<double> &less, vector<double> &more, double TID) { // how to get the data and plot it. Still to do: make histogram beatiful ;-)
	PixDBData data("scan", (fname + ":/" + sname + "/" + gname).c_str(), mname.c_str());

	TCanvas *c1=new TCanvas((fname + ":/" + sname + "/" + gname).c_str(),(fname + ":/" + sname + "/" + gname).c_str(),800,800);
	TH2F * myData = data.GetMap(0, PixScan::OCCUPANCY, -1);
	gStyle->SetOptStat(0);
	//cout << endl << "got data" << endl << endl;
	TH2F * myPlot = new TH2F(*(myData)/**((TH2F*)data.GetMap(0, PixScan::SCURVE_MEAN,0))*/);
	//myPlot->GetZaxis()->SetRangeUser(0,400);
	//cout << endl << "created TH2F" << (sname+".gif+5").c_str() << endl << endl;
	myPlot->SetName("");
	myPlot->SetTitle("");
	myPlot->Draw("COLZ");
	
	int numgood = 0, numless = 0, nummore = 0;
	TLatex *l1 = new TLatex();
	stringstream TID_string;
	if(TID > 0) TID_string << "TID "<<TID<<" Mrad"; 
	else TID_string << "TID "<<0<<" Mrad";
	l1->DrawText(.85, .55,TID_string.str().c_str());

	for (int numcol = 0; numcol < myPlot->GetNbinsX(); numcol++) {
		for (int numrow = 0; numrow < myPlot->GetNbinsY(); numrow++) {
			if ((myPlot->GetBinContent(numcol+1, numrow+1)) == 200) numgood++;
			else if ((myPlot->GetBinContent(numcol+1, numrow+1)) < 200)numless++;
			else if ((myPlot->GetBinContent(numcol+1, numrow+1)) > 200)nummore++;
		}
	}
	//c1->Print((sname+".gif+10").c_str()); //TO CREATE ANIMATED GIF - TAKE A LOT OF TIME
	c1->Delete();	
	good.push_back(numgood);
	less.push_back(numless);
	more.push_back(nummore);
	myPlot->Delete();
	//myData->Delete();
	//TApplication app("app", NULL, NULL);
	//myPlot->Draw("COLZ");
	//cout << "Press ctrl-C to stop" << endl;
	//app.Run();
}

void plotTOT(string fname, string sname, string gname, string mname, vector<double> &tot_mean, vector<double> &tot_rms){
	PixDBData data("scan", (fname + ":/" + sname + "/" + gname).c_str(), mname.c_str());
	TH2F * myData = data.GetMap(0, PixScan::TOT_MEAN, -1);
	TH2F * myPlot = new TH2F(*(myData)/**((TH2F*)data.GetMap(0, PixScan::SCURVE_MEAN,0))*/);
	TH1F *TOT_HISTO = new TH1F((fname + ":/" + sname + "/" + gname+"_HISTO").c_str(),(fname + ":/" + sname + "/" + gname+"_HISTO").c_str(),45,-0.5,14.5);
	for (int numcol = 0; numcol < myPlot->GetNbinsX(); numcol++) {
		for (int numrow = 0; numrow < myPlot->GetNbinsY(); numrow++) {
			//cout<<myPlot->GetBinContent(numcol+1, numrow+1)<<endl;
			
			if(myPlot->GetBinContent(numcol+1, numrow+1)>0) TOT_HISTO->Fill(myPlot->GetBinContent(numcol+1, numrow+1));
		}
	}
	tot_mean.push_back(TOT_HISTO->GetMean());
	tot_rms.push_back(TOT_HISTO->GetRMS());
	//for (int k = 0; k < tot_mean.size(); k++) {cout <<tot_mean[k]<<" "<<tot_rms[k]<<endl;}

	TOT_HISTO->Delete();
}

int returnsecond(std::string date, struct tm *tm) {


	//cout<<date<<endl;

	int month = atoi(date.substr(0, 2).c_str());
	int day = atoi(date.substr(3, 2).c_str());
	int year = atoi(date.substr(6, 4).c_str());
	int hour = atoi(date.substr(11, 2).c_str());
	int minute = atoi(date.substr(14, 2).c_str());
	int second = atoi(date.substr(17, 2).c_str());

	tm->tm_mday	= day;
	tm->tm_mon	= month - 1;
	tm->tm_year	= year - 1900;
	tm->tm_hour	= hour;
	tm->tm_min	= minute;
	tm->tm_sec	= second;

	return 1;
}

void Title(std::string sname,std::string name_string){
	TCanvas *c1=new TCanvas(name_string.c_str(),name_string.c_str(),800,800);
	TLatex *l1 = new TLatex();
	l1->DrawText(.15, .55,name_string.c_str());

	//c1->Print((sname+".gif+200").c_str());
	c1->Delete();	


}

void fill_vector(std::string name, int num_files, vector<double> &v_i_d_uncfg, vector<double>  &v_i_d_cfg, vector<double> &v_i_a_cfg, vector<double>  &v_time, vector<double>  &v_i_tot_cfg, vector<double> &v_tens_d_uncfg, vector<double> &v_tens_d_cfg, vector<double> &v_tens_a_cfg, vector<double> &good_digital, vector<double> &less_digital, vector<double> &more_digital, vector<double> &good_analog, vector<double> &less_analog, vector<double> &more_analog, TGraph* Charge_Time,vector<double> &tot_mean,vector<double> &tot_rms) {
	double i_d_uncfg, i_d_cfg, i_a_cfg, tens_d_uncfg, tens_d_cfg, tens_a_cfg;
	double utc_time;

	struct tm *tm;

	time_t rawtime;
	time ( &rawtime );
	tm = localtime ( &rawtime );

	time_t epoch;

	//strptime("03/02/2016_12:27:43", "%d-%m-%Y_%H:%M:%S", &tm);
	epoch = mktime(tm);



	for (int num = 0; num < num_files; num++) {
		std::string fname = name + to_string(num + 1) + ".root";
		cout << fname << endl;

		//if(num>=80) fname = name+"two_"+ to_string(num-79)+".root";
		//if(num>=26) fname ="IBLB01_2irr_afterAnnealing_GOOD_"+ to_string(num-25)+".root";
		vector<double> v;
		vector<string> v_name;
		string my_time;
		listScans(fname.c_str());
		fillData(fname.c_str(), my_time , v, v_name);
		//plotFDACd(fname.c_str())
		//cout << "DCS File Name: " << fname << endl;
		//cout << "DCS Time: " << my_time << endl;
		//cout << "DCS NAME: 		" << v_name[0] << " $ " << v_name[1] << " $ " << v_name[2] << " $ " << v_name[3] << " $ " << v_name[4] << " $ " << v_name[5] << endl;
		//cout << "DCS READING:	" << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << " " << v[4] << " " << v[5] << endl;
		//outputfile << my_time << " " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << " " << v[4] << " " << v[5] << endl;
		//cout << ScanName << " " << GroupName << " " << ModuleName << endl;

		string timestamp = my_time.substr(12, 10) + "_" + my_time.substr(0, 8);
		//cout << "DCS Time2: " << timestamp << endl;

		//string timestamp = my_time.substr(0, 8);

		//strptime(timestamp.c_str(), "%H:%M:%S", &tm)
		//if ( strptime(timestamp.c_str(), "%H:%M:%S", &tm) != NULL ) {
		if (returnsecond(timestamp, tm)) {
			//strptime(timestamp.c_str(), "%d-%m-%Y_%H:%M:%S", &tm);

			epoch = mktime(tm);

			//v_i_d_uncfg.push_back(i_d_uncfg);
			//v_i_d_cfg.push_back(i_d_cfg);
			//v_i_a_cfg.push_back(i_a_cfg);
			//v_i_tot_cfg.push_back(i_a_cfg + i_d_cfg);

			// V_VDDD2 uncfg I_VDDD2 uncfg V_VDDD2 cfg V_VDDA2 cfg I_VDDD cfg I_VDDA cfg

			for (int name_int = 0; name_int < v_name.size(); name_int++) {
				if (v_name[name_int] == "V_VDDD2 uncfg") v_tens_d_uncfg.push_back(v[name_int]);
				if (v_name[name_int] == "I_VDDD2 uncfg") v_i_d_uncfg.push_back(v[name_int]);
				if (v_name[name_int] == "V_VDDD2 cfg") v_tens_d_cfg.push_back(v[name_int]);
				if (v_name[name_int] == "I_VDDD cfg") {
					v_i_d_cfg.push_back(v[name_int]);
					i_d_cfg = v[name_int];
				}
				if (v_name[name_int] == "V_VDDA2 cfg") v_tens_a_cfg.push_back(v[name_int]);
				if (v_name[name_int] == "I_VDDA cfg") {
					v_i_a_cfg.push_back(v[name_int]);
					i_a_cfg = v[name_int];
				}
			}
			v_i_tot_cfg.push_back(i_a_cfg + i_d_cfg);

			if (v_time.size() > 0 && epoch < v_time[v_time.size() - 1] ) { cout << "CACCCCAAAAAA   " << epoch << " " << v_time[v_time.size() - 1] << endl; }
			//cout << "DCS TimeSec  : " << tm->tm_hour << " " << tm->tm_min << " " << tm->tm_sec << " " << endl;


			v_time.push_back((double)epoch);
			//cout<< setprecision(20) <<"FE "<< epoch <<endl;

		}

		//cout << "DCS READING SIZE:	" << v_i_d_uncfg.size() << " " << v_tens_d_uncfg.size() << " " << v_tens_d_cfg.size() << " " << v_i_d_cfg.size() << " " << v_tens_a_cfg.size() << " " << v_i_a_cfg.size() << endl;

		double TID = TIDValue(Charge_Time->Eval((double)epoch));
		cout<<TID<<endl;
		plotANALOGDIGITAL(fname.c_str(), "DIGITAL_TEST", GroupName, ModuleName, good_digital, less_digital, more_digital, TID);
		plotANALOGDIGITAL(fname.c_str(), "ANALOG_TEST", GroupName, ModuleName, good_analog, less_analog, more_analog, TID);
		plotTOT(fname.c_str(), "TOT_VERIF", GroupName, ModuleName,tot_mean,tot_rms);
		v.clear();
		v_name.clear();

	}

}

void fill_beamdump(std::string name, vector<double> &vx, vector<double> &vy, vector<double> &vz, double min, double max) {


	// myReadFile.open("2016_3_2__CurrentMeasWithChip.txt");
	// myReadFile_chip_current.open("First_irr.txt");
	ifstream myReadFile;

	struct tm *tm;
	time_t epoch;

	time_t rawtime;
	time ( &rawtime );
	tm = localtime ( &rawtime );


	myReadFile.open(name.c_str());
	double value;
	string date, hours;

	if (myReadFile.is_open()) {
		int counter = 0;
		while (!myReadFile.eof()) {
			myReadFile >> date >> hours >> value;
			//cout<< date << " " << hours << " " << value << endl;

			string timestamp = date + "_" + hours.substr(0, 8);

			double utc_time;
			int millisec = atoi(hours.substr(9, 3).c_str());
			// cout<<millisec/1000.<<endl;


			if ( returnsecond(timestamp, tm)) {
				//epoch = mktime(&tm);
				counter++;
				epoch = mktime(tm);

				utc_time = (double)epoch + (double)(millisec / 1000.);
				//cout<<setprecision(20)<<" "<< tm->tm_hour<<" "<<utc_time<<endl;

				if ((value > min && value < max) && utc_time > 1454495263.9509999752) {
					// if(true){
					if (vy.size()) {
						vz.push_back(vz[vz.size() - 1] + ((value + vy[vy.size() - 1]) * (utc_time - (vx[vx.size() - 1])) / 2 ));
					}
					else {
						vz.push_back(0);
					}
					vy.push_back(value);
					if (vx.size() > 0 && utc_time < vx[vx.size() - 1] ) cout << "CACCCCAAAAAA   " << utc_time << " " << vx[vx.size() - 1] << endl;
					vx.push_back(utc_time);

				}
				//  std::cout<<selection.str()<<std::endl;
				//std::cout<<timestamp<<"\t"<<epoch<<std::endl;
				// }
			}


			// if (value != 0)histo1->Fill(value);
		}
	}


}




void Fill_Charge_TID(vector<double> &v_time, vector<double> &v_Charge, vector<double> &v_TID, vector<double> &vx, TGraph *Charge_Time ,	int &index) {

	for (int i = 0; i < v_time.size(); i++) {

		if (v_time[i] > vx[0]) {
			v_Charge.push_back(Charge_Time->Eval(v_time[i]));
			v_TID.push_back((TIDValue(Charge_Time->Eval(v_time[i]))));
		}
		else index = i + 1;
	}
	
	cout << index << endl;

}

vector<double> offset(vector<double> a) {
	std::vector<double> v;
	double b = a[0]  ;

	for (int i = 0; i < a.size(); i++) {
		if (a[i] - b < 0) b = a[i];
		v.push_back(a[i] - b);

	}
	return v;
}

int main(int argc, char **argv) {
	//std::string name = "Data/irradiation_test_cyclotrone_";
	//std::string name = "Data/IBLB01_second_irradiation/IBLB01_2irr_beforeAnnealing_GOOD_";
	//std::string name = "Data/IBLB01_second_irradiation/IBLB01_2irr_afterAnnealing_GOOD_";

	vector<double> v_i_d_uncfg, v_i_d_cfg, v_i_a_cfg, v_time, v_i_tot_cfg;
	vector<double> v_tens_d_uncfg, v_tens_d_cfg, v_tens_a_cfg;
	vector<double> good_digital , less_digital , more_digital;
	vector<double> good_analog , less_analog , more_analog;	
	vector<double> tot_mean,tot_rms;

	vector<double> vx;
	vector<double> vy;
	vector<double> vz;

	cout<<"ciao"<<endl;
	fill_beamdump("Data/2016_3_2__CurrentMeasWithChip.txt", vx, vy, vz, -1e-12, 3e-9);
	fill_beamdump("Data/IBL01B_beam_dump_20160309.txt", vx, vy, vz, -1e-12, 3e-9);


	TGraph *Current_Time = new TGraph(vx.size(), &(vx[0]), &(vy[0]));
	TGraph *Charge_Time = new TGraph(vz.size(), &(vx[0]), &(vz[0]));



	std::string name = "Data/irradiation_test_cyclotrone_";
	fill_vector(name, 80, v_i_d_uncfg, v_i_d_cfg, v_i_a_cfg, v_time, v_i_tot_cfg, v_tens_d_uncfg, v_tens_d_cfg, v_tens_a_cfg , good_digital , less_digital , more_digital,good_analog , less_analog , more_analog, Charge_Time,tot_mean,tot_rms);
	name = "Data/irradiation_test_cyclotrone_two_";
	fill_vector(name,80,v_i_d_uncfg, v_i_d_cfg, v_i_a_cfg, v_time, v_i_tot_cfg,v_tens_d_uncfg, v_tens_d_cfg, v_tens_a_cfg , good_digital , less_digital , more_digital,good_analog , less_analog , more_analog, Charge_Time,tot_mean,tot_rms);


	//for (int k = 0; k < good_digital.size() ; k++) cout << good_digital[k] << " " << less_digital[k] << " " << more_digital[k] << " " << good_digital[k] + less_digital[k] + more_digital[k] << endl;

	//Title("ANALOG_TEST","One week annealing without BEAM");
	//Title("DIGITAL_TEST","One week annealing without BEAM");
	name ="Data/IBLB01_second_irradiation/IBLB01_2irr_beforeAnnealing_GOOD_";
	//25
	fill_vector(name, 25, v_i_d_uncfg, v_i_d_cfg, v_i_a_cfg, v_time, v_i_tot_cfg, v_tens_d_uncfg, v_tens_d_cfg, v_tens_a_cfg , good_digital , less_digital , more_digital,good_analog , less_analog , more_analog, Charge_Time,tot_mean,tot_rms);
	//Title("DIGITAL_TEST","FE OFF - BEAM ON");
	//Title("ANALOG_TEST","FE OFF - BEAM ON");
	name = "Data/IBLB01_second_irradiation/IBLB01_2irr_afterAnnealing_GOOD_";
	fill_vector(name, 22, v_i_d_uncfg, v_i_d_cfg, v_i_a_cfg, v_time, v_i_tot_cfg, v_tens_d_uncfg, v_tens_d_cfg, v_tens_a_cfg , good_digital , less_digital , more_digital,good_analog , less_analog , more_analog, Charge_Time,tot_mean,tot_rms);
	name = "Data/IBLB01_second_irradiation/IBLB01_2irr_afterAnnealing_BEAM_OFF_";
	fill_vector(name, 19, v_i_d_uncfg, v_i_d_cfg, v_i_a_cfg, v_time, v_i_tot_cfg, v_tens_d_uncfg, v_tens_d_cfg, v_tens_a_cfg , good_digital , less_digital , more_digital,good_analog , less_analog , more_analog, Charge_Time,tot_mean,tot_rms);

	//for (int k = 0; k < v_i_d_uncfg.size(); k++) {cout <<setprecision(20)<<v_time[k]<< " DCS READING:	" << v_i_d_uncfg[k] << " " << v_tens_d_uncfg[k] << " " << v_tens_d_cfg[k] << " " << v_i_d_cfg[k] << " " << v_tens_a_cfg[k] << " " << v_i_a_cfg[k] << endl;}



	cout << setprecision(10) << vx[0] << " " << vy[0] << " " << vz[0] << endl;
	cout << setprecision(10) << vx[1] << " " << vy[1] << " " << vz[1] << endl;


	TGraph *Good_Digital_Time = new TGraph(v_time.size(), &(v_time[0]), &(good_digital[0]));
	TGraph *Less_Digital_Time = new TGraph(v_time.size(), &(v_time[0]), &(less_digital[0]));
	TGraph *More_Digital_Time = new TGraph(v_time.size(), &(v_time[0]), &(more_digital[0]));

	Good_Digital_Time->SetName("Good_Digital_Time");
	Less_Digital_Time->SetName("Less_Digital_Time");
	More_Digital_Time->SetName("More_Digital_Time");

	TGraph *Good_Analog_Time = new TGraph(v_time.size(), &(v_time[0]), &(good_analog[0]));
	TGraph *Less_Analog_Time = new TGraph(v_time.size(), &(v_time[0]), &(less_analog[0]));
	TGraph *More_Analog_Time = new TGraph(v_time.size(), &(v_time[0]), &(more_analog[0]));

	Good_Analog_Time->SetName("Good_Analog_Time");
	Less_Analog_Time->SetName("Less_Analog_Time");
	More_Analog_Time->SetName("More_Analog_Time");

	TGraphErrors *TOTMEAN_Time = new TGraphErrors(v_time.size(), &(v_time[0]), &(tot_mean[0]),0,&(tot_rms[0]));
	//TGraph *TOTRMS_Time  = new TGraph(v_time.size(), &(v_time[0]), &(tot_rms[0]));

	TGraph *Iddd_uncfg_Time = new TGraph(v_time.size(), &(v_time[0]), &(v_i_d_uncfg[0]));
	TGraph *Iddd_cfg_Time = new TGraph(v_time.size(), &(v_time[0]), &(v_i_d_cfg[0]));
	TGraph *Idda_cfg_Time = new TGraph(v_time.size(), &(v_time[0]), &(v_i_a_cfg[0]));
	TGraph *Itot_cfg_Time = new TGraph(v_time.size(), &(v_time[0]), &(v_i_tot_cfg[0]));

	Iddd_cfg_Time->SetTitle("Iddd_cfg_Time");
	Idda_cfg_Time->SetTitle("Idda_cfg_Time");
	Iddd_uncfg_Time->SetTitle("Iddd_uncfg_Time");
	Itot_cfg_Time->SetTitle("Itot_cfg_Time");

	TOTMEAN_Time->SetTitle("TOTMEAN_Time");
	TOTMEAN_Time->SetName("TOTMEAN_Time");


	Iddd_cfg_Time->GetYaxis()->SetTitle("FE-I4 Current [A]");
	Idda_cfg_Time->GetYaxis()->SetTitle("FE-I4 Current [A]");
	Iddd_uncfg_Time->GetYaxis()->SetTitle("FE-I4 Current [A]");
	Itot_cfg_Time->GetYaxis()->SetTitle("FE-I4 Current [A]");

	Iddd_cfg_Time->SetName("Iddd_cfg_Time");
	Idda_cfg_Time->SetName("Idda_cfg_Time");
	Iddd_uncfg_Time->SetName("Iddd_uncfg_Time");
	Itot_cfg_Time->SetName("Itot_cfg_Time");


	vector<double> v_i_d_uncfg_offset, v_i_d_cfg_offset, v_i_a_cfg_offset, v_i_tot_cfg_offset;

	v_i_d_uncfg_offset = offset(v_i_d_uncfg);
	v_i_d_cfg_offset = offset(v_i_d_cfg);
	v_i_a_cfg_offset = offset(v_i_a_cfg);
	v_i_tot_cfg_offset = offset(v_i_tot_cfg);

	TGraph *Iddd_uncfg_Time_offset = new TGraph(v_time.size(), &(v_time[0]), &(v_i_d_uncfg_offset[0]));
	TGraph *Iddd_cfg_Time_offset = new TGraph(v_time.size(), &(v_time[0]), &(v_i_d_cfg_offset[0]));
	TGraph *Idda_cfg_Time_offset = new TGraph(v_time.size(), &(v_time[0]), &(v_i_a_cfg_offset[0]));
	TGraph *Itot_cfg_Time_offset = new TGraph(v_time.size(), &(v_time[0]), &(v_i_tot_cfg_offset[0]));

	vector<double> v_Charge;
	vector<double> v_TID;
	int index = 0;

	Fill_Charge_TID(v_time, v_Charge, v_TID, vx, Charge_Time, index);
	// vector<double> v_i_d_uncfg_rescale, v_i_d_cfg_rescale, v_i_a_cfg_rescale, v_Charge_rescale, v_i_tot_cfg_rescale;
	cout << index << endl;

	TGraph *Iddd_uncfg_BeamCurrent = new TGraph(v_Charge.size(), &(v_Charge[0]), &(v_i_d_uncfg[index]));
	TGraph *Iddd_cfg_BeamCurrent = new TGraph(v_Charge.size(), &(v_Charge[0]), &(v_i_d_cfg[index]));
	TGraph *Idda_cfg_BeamCurrent = new TGraph(v_Charge.size(), &(v_Charge[0]), &(v_i_a_cfg[index]));
	TGraph *Itot_cfg_BeamCurrent = new TGraph(v_Charge.size(), &(v_Charge[0]), &(v_i_tot_cfg[index]));

	TGraph *Iddd_uncfg_BeamCurrent_offset = new TGraph(v_Charge.size(), &(v_Charge[0]), &(v_i_d_uncfg_offset[index]));
	TGraph *Iddd_cfg_BeamCurrent_offset = new TGraph(v_Charge.size(), &(v_Charge[0]), &(v_i_d_cfg_offset[index]));
	TGraph *Idda_cfg_BeamCurrent_offset = new TGraph(v_Charge.size(), &(v_Charge[0]), &(v_i_a_cfg_offset[index]));
	TGraph *Itot_cfg_BeamCurrent_offset = new TGraph(v_Charge.size(), &(v_Charge[0]), &(v_i_tot_cfg_offset[index]));

	TGraph *Iddd_uncfg_TID = new TGraph(v_TID.size(), &(v_TID[0]), &(v_i_d_uncfg[index]));
	TGraph *Iddd_cfg_TID = new TGraph(v_TID.size(), &(v_TID[0]), &(v_i_d_cfg[index]));
	TGraph *Idda_cfg_TID = new TGraph(v_TID.size(), &(v_TID[0]), &(v_i_a_cfg[index]));
	TGraph *Itot_cfg_TID = new TGraph(v_TID.size(), &(v_TID[0]), &(v_i_tot_cfg[index]));

	TGraph *Iddd_uncfg_TID_offset = new TGraph(v_TID.size(), &(v_TID[0]), &(v_i_d_uncfg_offset[index]));
	TGraph *Iddd_cfg_TID_offset = new TGraph(v_TID.size(), &(v_TID[0]), &(v_i_d_cfg_offset[index]));
	TGraph *Idda_cfg_TID_offset = new TGraph(v_TID.size(), &(v_TID[0]), &(v_i_a_cfg_offset[index]));
	TGraph *Itot_cfg_TID_offset = new TGraph(v_TID.size(), &(v_TID[0]), &(v_i_tot_cfg_offset[index]));

	TGraph *Good_Digital_TID = new TGraph(v_TID.size(), &(v_TID[0]), &(good_digital[index]));
	TGraph *Less_Digital_TID = new TGraph(v_TID.size(), &(v_TID[0]), &(less_digital[index]));
	TGraph *More_Digital_TID = new TGraph(v_TID.size(), &(v_TID[0]), &(more_digital[index]));

	Good_Digital_TID->SetName("Good_Digital_TID");
	Less_Digital_TID->SetName("Less_Digital_TID");
	More_Digital_TID->SetName("More_Digital_TID");

	TGraph *Good_Analog_TID = new TGraph(v_TID.size(), &(v_TID[0]), &(good_analog[index]));
	TGraph *Less_Analog_TID = new TGraph(v_TID.size(), &(v_TID[0]), &(less_analog[index]));
	TGraph *More_Analog_TID = new TGraph(v_TID.size(), &(v_TID[0]), &(more_analog[index]));

	TGraphErrors *TOTMEAN_TID = new TGraphErrors(v_TID.size(), &(v_TID[0]), &(tot_mean[0]),0,&(tot_rms[0]));
	//TGraph *TOTRMS_TID  = new TGraph(v_TID.size(), &(v_TID[0]), &(tot_rms[0]));

	TOTMEAN_TID->SetTitle("TOTMEAN_TID");
	TOTMEAN_TID->SetName("TOTMEAN_TID");

	Good_Analog_TID->SetName("Good_Analog_TID");
	Less_Analog_TID->SetName("Less_Analog_TID");
	More_Analog_TID->SetName("More_Analog_TID");

	Iddd_cfg_Time_offset->SetName("Iddd_cfg_Time_offset");
	Idda_cfg_Time_offset->SetName("Idda_cfg_Time_offset");
	Iddd_uncfg_Time_offset->SetName("Iddd_uncfg_Time_offset");
	Itot_cfg_Time_offset->SetName("Itot_cfg_Time_offset");
	Current_Time->SetName("Current_Time");
	Charge_Time->SetName("Charge_Time");

	Iddd_cfg_BeamCurrent_offset->SetName("Iddd_cfg_BeamCurrent_offset");
	Idda_cfg_BeamCurrent_offset->SetName("Idda_cfg_BeamCurrent_offset");
	Iddd_uncfg_BeamCurrent_offset->SetName("Iddd_uncfg_BeamCurrent_offset");
	Itot_cfg_BeamCurrent_offset->SetName("Itot_cfg_BeamCurrent_offset");

	Iddd_cfg_BeamCurrent->SetName("Iddd_cfg_BeamCurrent");
	Idda_cfg_BeamCurrent->SetName("Idda_cfg_BeamCurrent");
	Iddd_uncfg_BeamCurrent->SetName("Iddd_uncfg_BeamCurrent");
	Itot_cfg_BeamCurrent->SetName("Itot_cfg_BeamCurrent");

	Iddd_cfg_TID_offset->SetName("Iddd_cfg_TID_offset");
	Idda_cfg_TID_offset->SetName("Idda_cfg_TID_offset");
	Iddd_uncfg_TID_offset->SetName("Iddd_uncfg_TID_offset");
	Itot_cfg_TID_offset->SetName("Itot_cfg_TID_offset");


	Iddd_cfg_TID->SetName("Iddd_cfg_TID");
	Idda_cfg_TID->SetName("Idda_cfg_TID");
	Iddd_uncfg_TID->SetName("Iddd_uncfg_TID");
	Itot_cfg_TID->SetName("Itot_cfg_TID");


	TFile *output_file = new TFile("FEI4_Irr_Bern_Cyclotron_10Mrad_Chip1.root", "RECREATE");

	Iddd_cfg_Time_offset->Write();
	Idda_cfg_Time_offset->Write();
	Iddd_uncfg_Time_offset->Write();
	Itot_cfg_Time_offset->Write();

	Current_Time->Write();
	Charge_Time->Write();

	TOTMEAN_Time->Write();
	//TOTRMS_Time->Write();
	TOTMEAN_TID->Write();
	//TOTRMS_TID->Write();
	
	Iddd_cfg_Time->Write();
	Idda_cfg_Time->Write();
	Iddd_uncfg_Time->Write();
	Itot_cfg_Time->Write();

	Iddd_cfg_BeamCurrent->Write();
	Idda_cfg_BeamCurrent->Write();
	Iddd_uncfg_BeamCurrent->Write();
	Itot_cfg_BeamCurrent->Write();

	Iddd_cfg_BeamCurrent_offset->Write();
	Idda_cfg_BeamCurrent_offset->Write();
	Iddd_uncfg_BeamCurrent_offset->Write();
	Itot_cfg_BeamCurrent_offset->Write();

	Iddd_cfg_TID->Write();
	Idda_cfg_TID->Write();
	Iddd_uncfg_TID->Write();
	Itot_cfg_TID->Write();

	Iddd_cfg_TID_offset->Write();
	Idda_cfg_TID_offset->Write();
	Iddd_uncfg_TID_offset->Write();
	Itot_cfg_TID_offset->Write();

	Good_Digital_Time->Write();
	Less_Digital_Time->Write();
	More_Digital_Time->Write();

	Good_Analog_Time->Write();
	Less_Analog_Time->Write();
	More_Analog_Time->Write();

	Good_Digital_TID->Write();
	Less_Digital_TID->Write();
	More_Digital_TID->Write();

	Good_Analog_TID->Write();
	Less_Analog_TID->Write();
	More_Analog_TID->Write();

	return 0;
}

