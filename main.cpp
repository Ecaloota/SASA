/* Parse an input .xyz file and return a .cube file containing
solvent-excluded surface area data encoded within a 
binary grid, for visualisation within IQMol. */

#include <fstream>
#include <iostream>
#include <tuple>
#include "boost/multi_array.hpp"
#include <QVector>
#include <QString>
#include <QStringList>
#include <QFile>
#include <QElapsedTimer>
#include <QTextStream>

//#include "QGLViewer/vec.h"

using namespace std;

typedef boost::multi_array<double, 3> Array3D;
typedef std::tuple<unsigned, unsigned, unsigned> GridPoint;

double const BohrRadius          = 5.2917721092e-11;
double const BohrToAngstrom      = BohrRadius*1.0e10;
double const AngstromToBohr      = 1.0/BohrToAngstrom;

// Need to work on integrating information about atomtypes
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

class Solvent {
	private:
		string m_solvent_name;
		double m_radius;

	public:
		Solvent() : m_solvent_name(), m_radius() {}

		Solvent(string solvent_name, double radius) : m_solvent_name(solvent_name), m_radius(radius) {}

		string getSolventName() {return m_solvent_name;}
		double getSolventRadius() {return m_radius;}
};

// Could work on improving speed
QVector<double> getBBMax(QVector<Atom> atomVector, Solvent solvent_instance) {
	
	float maxX = 0;
	float maxY = 0;
	float maxZ = 0;	
	float maxVDW = 0;
	double solvent_radius = solvent_instance.getSolventRadius();

	for(int iatom=0; iatom < atomVector.size(); iatom++) {

		if(atomVector[iatom].getvdW_radius(atomVector[iatom]) > maxVDW) { maxVDW = atomVector[iatom].getvdW_radius(atomVector[iatom]); }

		if(atomVector[iatom].getX() > maxX) {maxX = atomVector[iatom].getX();}
		if(atomVector[iatom].getY() > maxY) {maxY = atomVector[iatom].getY();}
		if(atomVector[iatom].getZ() > maxZ) {maxZ = atomVector[iatom].getZ();}
	}
	
	maxX += (maxVDW + solvent_radius); // largest VDW radius in input + some data points
	maxY += (maxVDW + solvent_radius); // fudge factor
	maxZ += (maxVDW + solvent_radius); // fudge factor

	QVector<double> bbmax = {maxX, maxY, maxZ};

	return bbmax;
}; 

// Could work on improving speed
QVector<double> getBBMin(QVector<Atom> atomVector, Solvent solvent_instance) {

	float minX = 0;
	float minY = 0;
	float minZ = 0;
	float maxVDW = 0;
	double solvent_radius = solvent_instance.getSolventRadius();

	for(int iatom=0; iatom < atomVector.size(); iatom++) {

		if(atomVector[iatom].getvdW_radius(atomVector[iatom]) > maxVDW) { maxVDW = atomVector[iatom].getvdW_radius(atomVector[iatom]); }

		if(atomVector[iatom].getX() < minX){minX = atomVector[iatom].getX();}
		if(atomVector[iatom].getY() < minY){minY = atomVector[iatom].getY();}
		if(atomVector[iatom].getZ() < minZ){minZ = atomVector[iatom].getZ();}
	}

	minX -= (maxVDW + solvent_radius); // largest VDW radius in input + some data points
	minY -= (maxVDW + solvent_radius); // fudge factor
	minZ -= (maxVDW + solvent_radius); // fudge factor

	QVector<double> bbmin = {minX, minY, minZ};

	return bbmin;
};

// Can make these refer to memory allocations rather than instances
float DistanceBetween3DPoints(float x1, float x2, float y1, float y2, float z1, float z2) {
	float distance = sqrt( pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2) );
	return distance;
}

QVector<Atom> parseXYZ( const string& infile ) {
	unsigned natoms;
	int atomIndex = 1;
	float xcoord, ycoord, zcoord;
	string atomtype; 
	QVector<Atom> atomVector;

	ifstream input_file(infile.c_str(), std::ifstream::in);
	if(!input_file.good()) {throw "The file does not exist! ";}

	input_file >> natoms;
	for(unsigned iatom=0; iatom<natoms; iatom++) {
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

GridData doSESTest(GridSize mygridsize, QVector<Atom> atomVector, GridData mygriddata, Solvent solvent) {

	/// Find the admissible points ///
	/// Rather, eliminate inadmissible points ///

	// GridData entries are initialised to 0, which is 'inadmissible'
	double solvent_radius = solvent.getSolventRadius();
	int maxNPWS = solvent_radius/mygridsize.delta()[0];
	QVector<GridPoint> admissiblePointVector;

	for (unsigned i=0; i<mygridsize.nx(); i++) {
		for (unsigned j=0; j<mygridsize.ny(); j++) {
			for (unsigned k=0; k<mygridsize.nz(); k++) {

				float x1 = mygridsize.origin()[0] + i*mygridsize.delta()[0];
				float y1 = mygridsize.origin()[1] + j*mygridsize.delta()[1];
				float z1 = mygridsize.origin()[2] + k*mygridsize.delta()[2];

				// Pretend point is admissible, then invalidate this if it is not true.
				mygriddata.operator()(i,j,k) = 1;

				for (int iatom=0; iatom<atomVector.size(); iatom++) {

					float x2 = atomVector[iatom].getX();
					float y2 = atomVector[iatom].getY();
					float z2 = atomVector[iatom].getZ();

					if (DistanceBetween3DPoints(x1, x2, y1, y2, z1, z2) < (atomVector[iatom].getvdW_radius(atomVector[iatom]) + solvent_radius)) {
						mygriddata.operator()(i,j,k) = 0;
						break;
					}
				}
			}
		}
	}

	// Store admissible points in memory
	// Alternatively, we could add all points in memory in first step and remove those
	// that end up failing the test, removing a second iteration through the whole grid
	// but increasing the number of memory allocations / de-allocations... ?
	for (unsigned i=0; i<mygridsize.nx(); i++) {
		for (unsigned j=0; j<mygridsize.ny(); j++) {
			for (unsigned k=0; k<mygridsize.nz(); k++) {

				if (mygriddata.operator()(i,j,k) == 1) {
					GridPoint admissiblePoint(i,j,k);
					admissiblePointVector.append(admissiblePoint);

				}
			}
		}
	}
	
	/// Find accessible points using the admissible subset ///
	// If the gridpoint refers to every atom and is outside each of them, do the rest of the test.
	// This involves assigning all points within the solvent sphere that is centred on (i,j,k) to 1

	for (int indx=0; indx<admissiblePointVector.size(); indx++) {

		unsigned i = get<0>(admissiblePointVector[indx]);
		unsigned j = get<1>(admissiblePointVector[indx]);
		unsigned k = get<2>(admissiblePointVector[indx]);

		float x1 = mygridsize.origin()[0] + i*mygridsize.delta()[0];
		float y1 = mygridsize.origin()[1] + j*mygridsize.delta()[1];
		float z1 = mygridsize.origin()[2] + k*mygridsize.delta()[2];

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
						mygriddata.operator()(ind1, ind2, ind3) = 1;
					}
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

	QElapsedTimer myTimer;
	myTimer.start();

	int Quality = 1;

	QVector<Atom> atomVector = parseXYZ("1ubq.xyz");
	QString filePath = "1ubq-SEStest-240817.cube";
	Solvent DefaultSolvent("Water", 1.4);

	QVector<double> BBMin = getBBMin(atomVector, DefaultSolvent);
	QVector<double> BBMax = getBBMax(atomVector, DefaultSolvent);

	GridSize myGridSize(BBMin, BBMax, Quality);
	GridData myGridData(myGridSize);	

	//GridData test = doVDWTest(myGridSize, atomVector, myGridData);
	GridData test = doSESTest(myGridSize, atomVector, myGridData, DefaultSolvent);

	QStringList coordinates = getCoordinates(atomVector);
	test.saveToCubeFile(filePath, coordinates, false);

	int testTimeNS = myTimer.elapsed(); // time taken in ns
	
	cout << " " << endl;
	cout << "SESTest completed in " << testTimeNS << " ns @ Quality = " << Quality << "." << endl;
	cout << " " << endl;

	return 0;
}




