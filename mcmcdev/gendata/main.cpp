#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

int main(){
	
	string outdatafile = "../testdata.dat";
	ofstream dataout;
	dataout.open(outdatafile);
	
	int ndps = 20;
	double m_fid = 1.0;
	double c_fid = 2.0;
	double m_off = 0.1;
	double c_off = 0.1;
	for(int dat = 0; dat < ndps; dat++){
		double x = dat * 0.1;
		double m = m_fid * ( 1.0 - m_off * rand()/(double)RAND_MAX );
		double c = c_fid * ( 1.0 - c_off * rand()/(double)RAND_MAX );
		double y = m * x + c;
		dataout << x << " " << y << endl;
	}
	
	dataout.close();
	
	
}