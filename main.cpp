/* Parse an input .xyz file and return a .cube file containing
solvent-accessible surface area data encoded within a 
binary array, for visualisation within IQMol. */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <stdio.h>

using namespace std;

struct Coordinates {
	float x_coord;
	float y_coord;
	float z_coord;
};

class Atom {
	private:
		string m_atomtype;
		Coordinates m_coordinates;
		int m_atom_index;
		float m_vdW_radius;
	public:
		
		Atom(string &atomtype, Coordinates coordinates, int atom_index, float vdW_radius) :
			m_atomtype(atomtype),m_coordinates(coordinates),m_atom_index(atom_index),m_vdW_radius(vdW_radius){}

		void printAtom() {
			cout << m_atom_index << setw(5) << m_atomtype << setw(7) << m_vdW_radius << setw(12) 
			<< m_coordinates.x_coord << setw(12) << m_coordinates.y_coord << setw(12) << m_coordinates.z_coord << endl;
		}

		float getX() {return m_coordinates.x_coord;}
		float getY() {return m_coordinates.y_coord;}
		float getZ() {return m_coordinates.z_coord;}

		float getvdW() {return m_vdW_radius;}
};

class Solvent {
	private:
		string m_solvent_name; // the name of the solvent (e.g. "h2o")
		float m_solvent_radius; // Or diameter, would save us a multiplication step down the track.
	public:
		Solvent(string solvent_name, float solvent_radius) : // This is our solvent Generator
			m_solvent_name(solvent_name),m_solvent_radius(solvent_radius) {
			}
		void printSolvent() {
			cout << "Solvent: " << m_solvent_name << " (" << m_solvent_radius << ")" << endl; // prints out: "m_solvent_name (m_solvent_radius)"
		}
		float getRadius() {return m_solvent_radius;}
};

class Griddata {
	private:
		vector<float> m_bbmin;
		vector<float>  m_bbmax;
		int m_quality;

	public:
		Griddata(vector<float> bbmin, vector<float> bbmax, int quality) : 
			m_bbmin(bbmin),m_bbmax(bbmax),m_quality(quality){}

		//int getNPX() {return m_npx;}
		//int getNPY() {return m_npy;}
		//int getNPZ() {return m_npz;}
};

vector<float> get_bbmax(vector<Atom> v_atoms) {
	
	float maximum_X = 0;
	float maximum_Y = 0;
	float maximum_Z = 0;	

	for(size_t iatom=0; iatom < v_atoms.size(); iatom++) {

		if(v_atoms[iatom].getX() > maximum_X) {maximum_X = v_atoms[iatom].getX();}
		if(v_atoms[iatom].getY() > maximum_Y) {maximum_Y = v_atoms[iatom].getY();}
		if(v_atoms[iatom].getZ() > maximum_Z) {maximum_Z = v_atoms[iatom].getZ();}
	}
	
	maximum_X += 2; // fudge factor
	maximum_Y += 2; // fudge factor
	maximum_Z += 2; // fudge factor

	float bbmax_v[3] = { maximum_X, maximum_Y, maximum_Z };
	vector<float> bbmax(&bbmax_v[0], &bbmax_v[0]+3);

	return bbmax;
}; 

vector<float> get_bbmin(vector<Atom> v_atoms) {

	float minimum_X = 0;
	float minimum_Y = 0;
	float minimum_Z = 0;

	for(size_t iatom=0; iatom < v_atoms.size(); iatom++) {
		
		if(v_atoms[iatom].getX() < minimum_X){minimum_X = v_atoms[iatom].getX();}
		if(v_atoms[iatom].getY() < minimum_Y){minimum_Y = v_atoms[iatom].getY();}
		if(v_atoms[iatom].getZ() < minimum_Z){minimum_Z = v_atoms[iatom].getZ();}
	}

	minimum_X -= 2; // fudge factor
	minimum_Y -= 2; // fudge factor
	minimum_Z -= 2; // fudge factor

	float bbmin_v[3] = { minimum_X, minimum_Y, minimum_Z };
	vector<float> bbmin(&bbmin_v[0], &bbmin_v[0]+3);

	return bbmin;
};

// Obtain the vdw radii of the atoms found in the input xyz, as part of parse_xyz
float getvdW_radius(string atomtype) {
	if(atomtype=="H") {return 1.20;}
	if(atomtype=="C") {return 1.70;}
	if(atomtype=="O") {return 1.52;}
	else {return 1;}
};

// At this point, this returns a vector of Atoms.
vector<Atom> parse_xyz( const string& infile ) 
{
	ifstream input_file(infile.c_str(), std::ifstream::in);
	if(!input_file.good()) {throw "The file does not exist! ";}

	size_t natoms; 
	string atomtype; 
	Coordinates coordinates; 
	int atom_index; 
	vector<Atom> v_atoms;

	input_file >> natoms;
	for(size_t iatom = 0; iatom < natoms; iatom++) {
		atom_index = iatom;
		input_file >> atomtype >> coordinates.x_coord >> coordinates.y_coord >> coordinates.z_coord;
		v_atoms.push_back(Atom(atomtype, coordinates, atom_index, getvdW_radius(atomtype)));
	};
	
	return v_atoms;
}

int main( ) 
{
	char input_xyz[50] = "methanol.xyz";
	vector<Atom> v_atoms;
	vector<float> bbmax;
	vector<float> bbmin;
	float quality = 0.01;

	//cout << "Please provide an input (xyz) file: "; 
//	cin >> input_xyz;

	// Not the best error handling I've seen, but it's something...
	try {v_atoms = parse_xyz(input_xyz);} catch (const char* msg) {cerr << msg << endl;}
	
	// This could of course be a single function.
	bbmax = get_bbmax(v_atoms);
	bbmin = get_bbmin(v_atoms);

	Griddata mygrid(bbmin,bbmax,quality);

	const int npx = (bbmax[0]-bbmin[0])/quality;
	const int npy = (bbmax[1]-bbmin[1])/quality;
	const int npz = (bbmax[2]-bbmin[2])/quality;

	// Make a vector of npx by npy by npz dimensions filled with 1s.
	vector<vector<vector<int> > > myfullvector (npx,vector<vector<int> >(npy,vector <int>(npz,1)));

	// We are only outputting zeroes....

	for (int i = 0; i < npx; i++) {
		for (int j = 0; j < npy; j++) {
			for (int k = 0; k < npz; k++) {
				float gridpoint[3] = { (bbmin[0] + (i*quality)), (bbmin[1] + (j*quality)), (bbmin[2] + (k*quality)) };
				//cout << "X = " << gridpoint[0] << " Y = " << gridpoint[1] << " Z = " << gridpoint[2] << endl;
				for(int iatom=0; iatom < v_atoms.size(); iatom++) {
					if(sqrt(pow((gridpoint[0] - v_atoms[iatom].getX()),2.0) + pow((gridpoint[1] - v_atoms[iatom].getY()),2.0) + pow((gridpoint[2] - v_atoms[iatom].getZ()),2.0)) < v_atoms[iatom].getvdW()) {
						//cout << "i = " << i << " j = " << j << " k = " << k << " VALUE " << sqrt(pow((gridpoint[0] - v_atoms[iatom].getX()),2.0) + pow((gridpoint[1] - v_atoms[iatom].getY()),2.0) + pow((gridpoint[2] - v_atoms[iatom].getZ()),2.0)) << " VDW " << v_atoms[iatom].getvdW() << endl;
						myfullvector[i][j][k] = 0;
						break;
					}
				}
			}
    	}
	}

	// Successfully prints out the atom attributes
	//for(size_t iatom=0; iatom < v_atoms.size(); iatom++) {
	//	v_atoms[iatom].printAtom();
	//}

	//float h2o_radius = 1.4;
	//string h2o_name = "Water";
	//Solvent h2o(h2o_name, h2o_radius);
	//h2o.printSolvent();

	cout << "NPX = " << npx << endl;
	cout << "NPY = " << npy << endl;
	cout << "NPZ = " << npz << endl;

	/*
	for(int i=0; i < npx; i++) {
		for(int j=0; j < npy; j++) {
			for(int k=0; k < npz; k++) {
				cout << myfullvector[i][j][k];
				if(k % 6 == 5)
					cout << "\n";
			}
			cout << "\n";
		}
	}
	*/

	//cout << "We are finally finished" << endl;

	return 0;
}























