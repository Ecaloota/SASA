/* Parse an input .xyz file and return a .cube file containing
solvent-accessible surface area data encoded within a 
binary array, for visualisation within IQMol. */

#include <fstream>
#include <iostream>
#include <tuple>
#include "boost/multi_array.hpp"
#include <QVector>
#include <QString>
#include <QStringList>
#include <QFile>
#include <QSet>
#include <QTextStream>

//#include "QGLViewer/vec.h"

using namespace std;

typedef boost::multi_array<double, 3> Array3D;
typedef std::tuple<int, int> Dyad;
typedef std::tuple<int, int, int> Tryad;

double const BohrRadius          = 5.2917721092e-11;
double const BohrToAngstrom      = BohrRadius*1.0e10;
double const AngstromToBohr      = 1.0/BohrToAngstrom;

// An atom contains an atomType, a set of coords and an index.
class Atom {
	private:
		unsigned m_index;
		string m_atomtype;
		float m_xcoord;
		float m_ycoord;
		float m_zcoord;
	public:
		
		Atom(unsigned index, string &atomtype, float xcoord, float ycoord, float zcoord) :
			m_index(index),m_atomtype(atomtype),m_xcoord(xcoord), m_ycoord(ycoord), m_zcoord(zcoord) {}

		Atom() : m_index(), m_atomtype(), m_xcoord(), m_ycoord(), m_zcoord() {}

		string getAtomType() {return m_atomtype;}
		int getAtomIndex() {return m_index;}
		float getX() {return m_xcoord;}
		float getY() {return m_ycoord;}
		float getZ() {return m_zcoord;}

		bool operator==(const Atom &other) const
		{
			return m_index == other.m_index;
		}

		//bool operator<(const Atom &other) const {
		//	if (m_xcoord < other.m_xcoord) {return true;}
		//	else if (m_xcoord == other.m_xcoord) {
		//		if (m_ycoord < other.m_ycoord) {return true;}
		//		else if (m_ycoord == other.m_ycoord) {
		//			if (m_zcoord < other.m_zcoord) {return true;}
		//		}
		//	}
		//	return false;
		//}

		// Refer AtomicProperty.C
		float getvdW_radius(Atom iatom) {
			if(iatom.getAtomType()=="H") {return 1.20;}
			if(iatom.getAtomType()=="He") {return 1.40;}
			if(iatom.getAtomType()=="Li") {return 1.82;}
			if(iatom.getAtomType()=="C") {return 1.70;}
			if(iatom.getAtomType()=="N") {return 1.55;}
			if(iatom.getAtomType()=="O") {return 1.52;}
			if(iatom.getAtomType()=="F") {return 1.47;}
			if(iatom.getAtomType()=="Ne") {return 1.54;}
			if(iatom.getAtomType()=="Na") {return 2.27;}
			if(iatom.getAtomType()=="Mg") {return 1.73;}
			if(iatom.getAtomType()=="Si") {return 2.10;}
			if(iatom.getAtomType()=="P") {return 1.80;}
			if(iatom.getAtomType()=="S") {return 1.80;}
			if(iatom.getAtomType()=="Cl") {return 1.75;}
			if(iatom.getAtomType()=="Ar") {return 1.88;}
			if(iatom.getAtomType()=="K") {return 2.75;}
			if(iatom.getAtomType()=="Ni") {return 1.63;}
			if(iatom.getAtomType()=="Cu") {return 1.40;}
			if(iatom.getAtomType()=="Zn") {return 1.39;}
			if(iatom.getAtomType()=="Ga") {return 1.87;}
			if(iatom.getAtomType()=="As") {return 1.85;}
			if(iatom.getAtomType()=="Se") {return 1.90;}
			if(iatom.getAtomType()=="Br") {return 1.85;}
			if(iatom.getAtomType()=="Kr") {return 2.02;}
			if(iatom.getAtomType()=="Pd") {return 1.63;}
			if(iatom.getAtomType()=="Ag") {return 1.72;}
			if(iatom.getAtomType()=="Cd") {return 1.58;}
			if(iatom.getAtomType()=="In") {return 1.93;}
			if(iatom.getAtomType()=="Sn") {return 2.17;}
			if(iatom.getAtomType()=="Te") {return 2.06;}
			if(iatom.getAtomType()=="I") {return 1.98;}
			if(iatom.getAtomType()=="Xe") {return 2.16;}
			if(iatom.getAtomType()=="Pt") {return 1.75;}
			if(iatom.getAtomType()=="Au") {return 1.66;}
			if(iatom.getAtomType()=="Tl") {return 1.96;}
			if(iatom.getAtomType()=="Pb") {return 2.02;}
			if(iatom.getAtomType()=="U") {return 1.86;}
			else {return 1;}
		};

		int getAtomicNumber(Atom iatom) {
			if(iatom.getAtomType()=="H") {return 1;}
			if(iatom.getAtomType()=="He") {return 2;}
			if(iatom.getAtomType()=="Li") {return 3;}
			if(iatom.getAtomType()=="Be") {return 4;}
			if(iatom.getAtomType()=="B") {return 5;}
			if(iatom.getAtomType()=="C") {return 6;}
			if(iatom.getAtomType()=="N") {return 7;}
			if(iatom.getAtomType()=="O") {return 8;}
			if(iatom.getAtomType()=="F") {return 9;}
			if(iatom.getAtomType()=="Ne") {return 10;}
			if(iatom.getAtomType()=="Na") {return 11;}
			if(iatom.getAtomType()=="Mg") {return 12;}
			if(iatom.getAtomType()=="Al") {return 13;}
			if(iatom.getAtomType()=="Si") {return 14;}
			if(iatom.getAtomType()=="P") {return 15;}
			if(iatom.getAtomType()=="S") {return 16;}
			if(iatom.getAtomType()=="Cl") {return 17;}
			if(iatom.getAtomType()=="Ar") {return 18;}
			if(iatom.getAtomType()=="K") {return 19;}
			if(iatom.getAtomType()=="Ca") {return 20;}
			else {return 1;}
		}
};

typedef std::tuple<Atom, Atom> AtomDyad; // A pair of Atom instances.
typedef std::tuple<Atom, Atom, Atom> AtomTryad; // A triplet of Atom instances.

double stepSize(unsigned const quality) {
   // These spacings are chosen so that each step uses roughly four times as
   // many points as the previous one.
   double stepSize(0.0);
   switch (quality) {
      case 0:   stepSize = 1.000000;  break;
      case 1:   stepSize = 0.629961;  break;
      case 2:   stepSize = 0.396850;  break;
      case 3:   stepSize = 0.250000;  break;
      case 4:   stepSize = 0.157490;  break;
      case 5:   stepSize = 0.099213;  break;
      case 6:   stepSize = 0.062500;  break;
      case 7:   stepSize = 0.039373;  break;
      case 8:   stepSize = 0.024803;  break;
      case 9:   stepSize = 0.015625;  break;
      default:  stepSize = 0.250000;  break;
   }
   return stepSize;
}

class GridSize {
	private:
		QVector<double> m_origin;
        QVector<double> m_delta;
        unsigned m_nx;
        unsigned m_ny;
        unsigned m_nz;

	public:
		GridSize(QVector<double> const& origin, QVector<double> const& delta, unsigned const nx, unsigned const ny, unsigned const nz) : 
		m_origin(origin), m_delta(delta), m_nx(nx), m_ny(ny), m_nz(nz) {}

        GridSize(QVector<double> const& min, QVector<double> const& max, int const quality) 
        {
         	double d(stepSize(quality));

         	m_origin = min;
         	m_delta.fill(d,3);

         	QVector<double> delta = {(max[0]-min[0]), (max[1]-min[1]), (max[2]-min[2])};
   			delta = {delta[0]/d, delta[1]/d, delta[2]/d};

   			m_nx = std::ceil(delta[0]) + 1;
   			m_ny = std::ceil(delta[1]) + 1;
   			m_nz = std::ceil(delta[2]) + 1;
         }

        QVector<double> const& delta() const { return m_delta; }
        QVector<double> const& origin() const { return m_origin; }

        unsigned nx() const { return m_nx; }
        unsigned ny() const { return m_ny; }
        unsigned nz() const { return m_nz; }
};

class GridData {
	private:

		Array3D m_data;
		QVector<double> m_origin;
		QVector<double> m_delta;
		
	public:
		GridData(GridSize const& size) : m_origin(size.origin()), m_delta(size.delta())
		{
			Array3D::extent_gen extents;
			m_data.resize(extents[size.nx()][size.ny()][size.nz()]);
		}

		void getNumberOfPoints(unsigned& nx, unsigned& ny, unsigned& nz) const
		{
			nx = m_data.shape()[0];
			ny = m_data.shape()[1];
			nz = m_data.shape()[2];
		}

		QVector<double> const& origin() const { return m_origin; }
        QVector<double> const& delta() const {return m_delta; }

        // Read only
		double const& operator()(unsigned const i, unsigned const j, unsigned const k) const
         {
            return m_data[i][j][k];
         }

         // Read and write
         double& operator()(unsigned const i, unsigned const j, unsigned const k)
         {
            return m_data[i][j][k];
         }

		bool saveToCubeFile(QString const& filePath, QStringList const& coordinates, 
			bool const invertSign) const
		{
			QFile file(filePath);
			if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
				return 0;

			QTextStream out(&file);

			QStringList header;
			header << "Cube file";
			header << "Generated using IQmol";
   
			QVector<double> delta = {AngstromToBohr*m_delta[0], AngstromToBohr*m_delta[1], AngstromToBohr*m_delta[2]};
   			QVector<double> origin = {AngstromToBohr*m_origin[0], AngstromToBohr*m_origin[1], AngstromToBohr*m_origin[2]};

   			unsigned nx, ny, nz;
   			getNumberOfPoints(nx, ny, nz);

   			header << QString("%1 %2 %3 %4").arg(coordinates.size(), 5)
                                   			.arg(origin[0], 13, 'f', 6)
                                   			.arg(origin[1], 13, 'f', 6)
                                   			.arg(origin[2], 13, 'f', 6); 

   			header << QString("%1 %2 %3 %4").arg(nx, 5)
                                   			.arg(delta[0], 13, 'f', 6)
                                   			.arg(0.0, 13, 'f', 6)
                                   			.arg(0.0, 13, 'f', 6);

   			header << QString("%1 %2 %3 %4").arg(ny, 5)
                                   			.arg(0.0, 13, 'f', 6)
                                   			.arg(delta[1], 13, 'f', 6)
                                   			.arg(0.0, 13, 'f', 6);

   			header << QString("%1 %2 %3 %4").arg(nz, 5)
                                   			.arg(0.0, 13, 'f', 6)
                                   			.arg(0.0, 13, 'f', 6)
                                   			.arg(delta[2], 13, 'f', 6);

   			header << coordinates;

   			QByteArray buffer;
   			buffer.append(header.join("\n"));
   			buffer.append("\n");
   			
   			out << buffer;
   			buffer.clear();

   			double w;
   			unsigned col(0);

   			for (unsigned i = 0; i < nx; ++i) {
       			for (unsigned j = 0; j < ny; ++j) {
           			for (unsigned k = 0; k < nz; ++k, ++col) {
               			w = m_data[i][j][k];
               			if (invertSign) w = -w; 
               			if (w >= 0.0) buffer += " ";
               			buffer += QString::number(w, 'E', 5); 
               			if (col == 5) {
                  			col = -1; 
                  			buffer += "\n";
               			}else {
                  			buffer += " ";
               			}   
               		}

               		out << buffer;
               		buffer.clear();
               	}
            }   

            buffer += "\n";
            out << buffer;
            file.flush();
            file.close();

            return true;
        }
};

QVector<double> getBBMax(QVector<Atom> atomVector, int quality) {
	
	float maxX = 0;
	float maxY = 0;
	float maxZ = 0;	
	float maxVDW = 0;

	for(int iatom=0; iatom < atomVector.size(); iatom++) {

		if(atomVector[iatom].getvdW_radius(atomVector[iatom]) > maxVDW) { maxVDW = atomVector[iatom].getvdW_radius(atomVector[iatom]); }

		if(atomVector[iatom].getX() > maxX) {maxX = atomVector[iatom].getX();}
		if(atomVector[iatom].getY() > maxY) {maxY = atomVector[iatom].getY();}
		if(atomVector[iatom].getZ() > maxZ) {maxZ = atomVector[iatom].getZ();}
	}
	
	maxX += (maxVDW + 15*stepSize(quality)); // largest VDW radius in input + some data points
	maxY += (maxVDW + 15*stepSize(quality)); // fudge factor
	maxZ += (maxVDW + 15*stepSize(quality)); // fudge factor

	QVector<double> bbmax = {maxX, maxY, maxZ};

	return bbmax;
}; 

QVector<double> getBBMin(QVector<Atom> atomVector, int quality) {

	float minX = 0;
	float minY = 0;
	float minZ = 0;
	float maxVDW = 0;

	for(int iatom=0; iatom < atomVector.size(); iatom++) {

		if(atomVector[iatom].getvdW_radius(atomVector[iatom]) > maxVDW) { maxVDW = atomVector[iatom].getvdW_radius(atomVector[iatom]); }

		if(atomVector[iatom].getX() < minX){minX = atomVector[iatom].getX();}
		if(atomVector[iatom].getY() < minY){minY = atomVector[iatom].getY();}
		if(atomVector[iatom].getZ() < minZ){minZ = atomVector[iatom].getZ();}
	}

	minX -= (maxVDW + 15*stepSize(quality)); // largest VDW radius in input + some data points
	minY -= (maxVDW + 15*stepSize(quality)); // fudge factor
	minZ -= (maxVDW + 15*stepSize(quality)); // fudge factor

	QVector<double> bbmin = {minX, minY, minZ};

	return bbmin;
};

float DistanceBetween3DPoints(float x1, float x2, float y1, float y2, float z1, float z2) {

	float distance = sqrt( pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2) );
	return distance;
}

QList<AtomDyad> getUniquePairs(QVector<Atom> &atomVector) {

	QList<AtomDyad> AtomDyadPairList;

	for (int iatom = 0; iatom < atomVector.size(); iatom++) {
		for (int jatom = 0; jatom < atomVector.size(); jatom++) {
			
			if (iatom == jatom) { continue; }

			float x1 = atomVector[jatom].getX();
			float x2 = atomVector[iatom].getX();

			float y1 = atomVector[jatom].getY();
			float y2 = atomVector[iatom].getY();
			
			float z1 = atomVector[jatom].getZ();
			float z2 = atomVector[iatom].getZ();

			float Radii_Sum1 = atomVector[jatom].getvdW_radius(atomVector[jatom]) + atomVector[iatom].getvdW_radius(atomVector[iatom]);


			if (DistanceBetween3DPoints(x1, x2, y1, y2, z1, z2) <= Radii_Sum1) {
				AtomDyad MyDyad(atomVector[iatom], atomVector[jatom]);
				AtomDyad MyDuplicate(atomVector[jatom], atomVector[iatom]);

				if (!AtomDyadPairList.contains(MyDuplicate)) {
					AtomDyadPairList.append(MyDyad);
				}
			}
		}
	}

	//for(int pair=0; pair<AtomDyadPairList.size(); pair++) {
	//	cout << "Dyad: " << (get<0>(AtomDyadPairList[pair])).getAtomIndex() << " " << (get<1>(AtomDyadPairList[pair])).getAtomIndex() << " " << endl;
	//}

	return AtomDyadPairList;
}

QList<AtomTryad> getUniqueTriples(QVector<Atom> atomVector, QList<AtomDyad> AtomPairList) {

	QList<AtomTryad> AtomTriplesList;

	for (int iatom = 0; iatom < atomVector.size(); iatom ++) {
		for (int ipair = 0; ipair < AtomPairList.size(); ipair++ ) {
			// if the distance between atom 1 of the pair and iatom is less than their VDW radii AND the distance between atom 2 of the 
			// pair and iatom is less than their VDW radii AND the distance between atom 1 and 2 of the pair is less than their VDW radii, 
			// THEN add (iatom, jatom, katom) to the triples list. iatom is bonded to both jatom and katom.

			int jatom = get<0>(AtomPairList[ipair]).getAtomIndex() - 1; // Because we 'imported' the pair list, we need to convert the indices.
			int katom = get<1>(AtomPairList[ipair]).getAtomIndex() - 1;

			float Radii_Sum1 = atomVector[jatom].getvdW_radius(atomVector[jatom]) + atomVector[iatom].getvdW_radius(atomVector[iatom]);
			float Radii_Sum2 = atomVector[katom].getvdW_radius(atomVector[katom]) + atomVector[iatom].getvdW_radius(atomVector[iatom]);
			float Radii_Sum3 = atomVector[jatom].getvdW_radius(atomVector[jatom]) + atomVector[katom].getvdW_radius(atomVector[katom]);
			
			float x1 = atomVector[jatom].getX();
			float x2 = atomVector[iatom].getX();
			float x3 = atomVector[katom].getX();

			float y1 = atomVector[jatom].getY();
			float y2 = atomVector[iatom].getY();
			float y3 = atomVector[katom].getY();

			float z1 = atomVector[jatom].getZ();
			float z2 = atomVector[iatom].getZ();
			float z3 = atomVector[katom].getZ();

			// We want only i,j,k. Brevity is the soul of wit.
			AtomTryad Duplicate1(atomVector[iatom], atomVector[katom], atomVector[jatom]); // We don't want i,k,j
			AtomTryad Duplicate2(atomVector[jatom], atomVector[iatom], atomVector[katom]); // We don't want j,i,k
			AtomTryad Duplicate3(atomVector[jatom], atomVector[katom], atomVector[iatom]); // We don't want j,k,i
			AtomTryad Duplicate4(atomVector[katom], atomVector[iatom], atomVector[jatom]); // We don't want k,i,j
			AtomTryad Duplicate5(atomVector[katom], atomVector[jatom], atomVector[iatom]); // We don't want k,j,i

			// Points 1 are Jatoms, Points 2 are Iatoms, Points 3 are Katoms.
			if (DistanceBetween3DPoints(x1, x2, y1, y2, z1, z2) <= Radii_Sum1 && DistanceBetween3DPoints(x3, x2, y3, y2, z3, z2) <= Radii_Sum2 && DistanceBetween3DPoints(x1, x3, y1, y3, z1, z3) <= Radii_Sum3) {
				if (atomVector[iatom].getAtomIndex() == atomVector[jatom].getAtomIndex() || atomVector[iatom].getAtomIndex() == atomVector[katom].getAtomIndex()) {continue;} // Disregard duplicates
				if (AtomTriplesList.contains(Duplicate1) || AtomTriplesList.contains(Duplicate2) || AtomTriplesList.contains(Duplicate3) || AtomTriplesList.contains(Duplicate4) || AtomTriplesList.contains(Duplicate5)) {continue;} // Disregard duplicates.
				AtomTryad MyTryad(atomVector[iatom], atomVector[jatom], atomVector[katom]);
				AtomTriplesList.append(MyTryad);
			}
		}
	}

	// Print out the Tryads.
	//for (int tryad = 0; tryad < AtomTriplesList.size(); tryad++) {
	//	cout << "Tryad: " << get<0>(AtomTriplesList[tryad]).getAtomIndex() << " " << get<1>(AtomTriplesList[tryad]).getAtomIndex() << " " << get<2>(AtomTriplesList[tryad]).getAtomIndex() << " " << endl;
	//}

	return AtomTriplesList;
}

QVector<Atom> parseXYZ( const string& infile ) {
	ifstream input_file(infile.c_str(), std::ifstream::in);
	if(!input_file.good()) {throw "The file does not exist! ";}

	unsigned natoms;
	int atomIndex = 1;
	float xcoord, ycoord, zcoord;
	string atomtype; 
	QVector<Atom> atomVector;

	input_file >> natoms;
	for(unsigned iatom = 0; iatom < natoms; iatom++) {
		input_file >> atomtype >> xcoord >> ycoord >> zcoord;
		atomVector.append(Atom(atomIndex, atomtype, xcoord, ycoord, zcoord));
		atomIndex += 1;
	};
	
	return atomVector;
}

GridData doVDWTest(GridSize mygridsize, QVector<Atom> atomVector, GridData mygriddata) {

	for (unsigned i=0; i<mygridsize.nx(); i++) {
		for (unsigned j=0; j<mygridsize.ny(); j++) {
			for (unsigned k=0; k<mygridsize.nz(); k++) {
				
				float x1 = mygridsize.origin()[0] + i*mygridsize.delta()[0];
				float y1 = mygridsize.origin()[1] + j*mygridsize.delta()[1];
				float z1 = mygridsize.origin()[2] + k*mygridsize.delta()[2];

				for (int iatom=0; iatom<atomVector.size(); iatom++) {

					float x2 = atomVector[iatom].getX();
					float y2 = atomVector[iatom].getY();
					float z2 = atomVector[iatom].getZ();

					if (DistanceBetween3DPoints(x1, x2, y1, y2, z1, z2) <= atomVector[iatom].getvdW_radius(atomVector[iatom])) {
						mygriddata.operator()(i,j,k) = 1;
						break;
					}
				}
			}
		}
	}

	return mygriddata;
}

GridData doSESTest(GridSize mygridsize, QVector<Atom> atomVector, GridData mygriddata, QList<AtomDyad> AtomPairList, QList<AtomTryad> AtomTriplesList) {

	// Plan is as follows:
	// for (point pi) {
	// for (pair in AtomPairList) {
	//		if ||A-Pi|| >= Ra + Rs && ||B-Pi|| >= Rb + Rs:
	//			// Pi is admissible
	//			// Every Pj such that ||Pi - Pj|| <= Rs is solvent accessible: 0 -> 1
	//			// Default is 0, which is inaccessible
	//			// if a point is found to be accessible, point = 1.

	float solvent_radius = 1.4;
	int maxNPWS = solvent_radius/mygridsize.delta()[0]; // Assuming delta is the same for all dimensions

	for (unsigned i=0; i<mygridsize.nx(); i++) {
		for (unsigned j=0; j<mygridsize.ny(); j++) {
			for (unsigned k=0; k<mygridsize.nz(); k++) {

				float x1 = mygridsize.origin()[0] + i*mygridsize.delta()[0];
				float y1 = mygridsize.origin()[1] + j*mygridsize.delta()[1];
				float z1 = mygridsize.origin()[2] + k*mygridsize.delta()[2];

				for (int ipair=0; ipair < AtomPairList.size(); ipair++) {

					// X2, Y2, Z2 are the coords of Atom A, (element 0) in the pair
					float Radii_Sum1 = get<0>(AtomPairList[ipair]).getvdW_radius(get<0>(AtomPairList[ipair])) + solvent_radius;
					float x2 = get<0>(AtomPairList[ipair]).getX();
					float y2 = get<0>(AtomPairList[ipair]).getY();
					float z2 = get<0>(AtomPairList[ipair]).getZ();

					// X3, Y3, Z3 are the coords of Atom B in the pair
					float Radii_Sum2 = get<1>(AtomPairList[ipair]).getvdW_radius(get<1>(AtomPairList[ipair])) + solvent_radius;
					float x3 = get<1>(AtomPairList[ipair]).getX();
					float y3 = get<1>(AtomPairList[ipair]).getY();
					float z3 = get<1>(AtomPairList[ipair]).getZ();					

					if (DistanceBetween3DPoints(x1, x2, y1, y2, z1, z2) >= Radii_Sum1 && DistanceBetween3DPoints(x1, x3, y1, y3, z1, z3) >= Radii_Sum2) {
						// Pi (x1, y1, z1) is admissible. All admissible points are accessible, therefore (x1, y1, z1) -> 1.
						mygriddata.operator()(i,j,k) = 1;

						for (int d1 = -maxNPWS; d1 < maxNPWS+1; d1++) {
							for (int d2 = -maxNPWS; d2 < maxNPWS+1; d2++) {
								for (int d3 = -maxNPWS; d3 < maxNPWS+1; d3++) {

									int ind1 = i + d1;
									int ind2 = j + d2;
									int ind3 = k + d3;

									// Unneccesarily testing the starting gridpoint
									if (d1 == 0 && d2 == 0 && d3 ==0) {
										continue; 
									}

									// Trying to index outside the grid itself (less than minimum).
									if (ind1 < 0 || ind2 < 0 || ind3 < 0) {
										continue;
									}

									// Trying to index outside the grid itself (greater than maximum) - Remember zero indexing.
									if (ind1 >= mygridsize.nx() || ind2 >= mygridsize.ny() || ind3 >= mygridsize.nz()) {
										continue;
									}
									
									float xj1 = mygridsize.origin()[0] + (ind1)*mygridsize.delta()[0];
									float yj1 = mygridsize.origin()[1] + (ind2)*mygridsize.delta()[1];
									float zj1 = mygridsize.origin()[2] + (ind3)*mygridsize.delta()[2];

									if ( DistanceBetween3DPoints(x1, xj1, y1, yj1, z1, zj1) <= solvent_radius) {
										// The gridpoint is solvent accessible.
										mygriddata.operator()(ind1, ind2, ind3) = 1;
									}
								}
							}
						}
					}
					//cout << get<0>(AtomPairList[ipair]).getAtomIndex() << " " << get<1>(AtomPairList[ipair]).getAtomIndex() << endl;
				}
			}
		}
	}

	return mygriddata;
}

QStringList getCoordinates(QVector<Atom> atomVector) {

	QStringList coordinates;
	for (int iatom = 0; iatom < atomVector.size(); iatom++) {
		coordinates << QString("%1 %2 %3 %4 %5").arg(atomVector[iatom].getAtomicNumber(atomVector[iatom]), 5)
												.arg(0.0, 13, 'f', 6)
												.arg(atomVector[iatom].getX(), 13, 'f', 6)
												.arg(atomVector[iatom].getY(), 13, 'f', 6)
												.arg(atomVector[iatom].getZ(), 13, 'f', 6);
	}

	return coordinates;
}

int main( ) {

	int quality = 5;

	QVector<Atom> atomVector = parseXYZ("hydrogen.xyz");
	QString filePath = "hydrogen_SEStest.cube";

	QVector<double> bbmin = getBBMin(atomVector, quality);
	QVector<double> bbmax = getBBMax(atomVector, quality);
	GridSize mygridsize(bbmin, bbmax, quality);
	GridData mygriddata(mygridsize);
	
	QList<AtomDyad> AtomPairList = getUniquePairs(atomVector);
	QList<AtomTryad> AtomTriplesList = getUniqueTriples(atomVector, AtomPairList);

	//GridData test = doVDWTest(mygridsize, atomVector, mygriddata);
	GridData test = doSESTest(mygridsize, atomVector, mygriddata, AtomPairList, AtomTriplesList);
	QStringList coordinates = getCoordinates(atomVector);

	test.saveToCubeFile(filePath, coordinates, false);

	return 0;
}




