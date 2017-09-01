/* Parse an input .xyz file and return a .cube file containing
solvent-excluded surface area data encoded within a 
binary array, for visualisation within IQMol. */

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

//#include <openbabel>
#include "QGLViewer/vec.h"

using namespace std;

typedef boost::multi_array<double, 3> Array3D;
typedef std::tuple<unsigned, unsigned, unsigned> GridPoint;

double const BohrRadius          = 5.2917721092e-11;
double const BohrToAngstrom      = BohrRadius*1.0e10;
double const AngstromToBohr      = 1.0/BohrToAngstrom;

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

		string type() const {return m_atomtype;}
		int getAtomIndex() const {return m_index;}
		float getX() const {return m_xcoord;}
		float getY() const {return m_ycoord;}
		float getZ() const {return m_zcoord;}

		bool operator==(const Atom &other) const
		{
			return m_index == other.m_index;
		}

		// Refer AtomicProperty.C
		float vdw(Atom iatom) {
			if(iatom.type()=="H") {return 1.20;}
			if(iatom.type()=="He") {return 1.40;}
			if(iatom.type()=="Li") {return 1.82;}
			if(iatom.type()=="C") {return 1.70;}
			if(iatom.type()=="N") {return 1.55;}
			if(iatom.type()=="O") {return 1.52;}
			if(iatom.type()=="F") {return 1.47;}
			if(iatom.type()=="Ne") {return 1.54;}
			if(iatom.type()=="Na") {return 2.27;}
			if(iatom.type()=="Mg") {return 1.73;}
			if(iatom.type()=="Si") {return 2.10;}
			if(iatom.type()=="P") {return 1.80;}
			if(iatom.type()=="S") {return 1.80;}
			if(iatom.type()=="Cl") {return 1.75;}
			if(iatom.type()=="Ar") {return 1.88;}
			if(iatom.type()=="K") {return 2.75;}
			if(iatom.type()=="Ni") {return 1.63;}
			if(iatom.type()=="Cu") {return 1.40;}
			if(iatom.type()=="Zn") {return 1.39;}
			if(iatom.type()=="Ga") {return 1.87;}
			if(iatom.type()=="As") {return 1.85;}
			if(iatom.type()=="Se") {return 1.90;}
			if(iatom.type()=="Br") {return 1.85;}
			if(iatom.type()=="Kr") {return 2.02;}
			if(iatom.type()=="Pd") {return 1.63;}
			if(iatom.type()=="Ag") {return 1.72;}
			if(iatom.type()=="Cd") {return 1.58;}
			if(iatom.type()=="In") {return 1.93;}
			if(iatom.type()=="Sn") {return 2.17;}
			if(iatom.type()=="Te") {return 2.06;}
			if(iatom.type()=="I") {return 1.98;}
			if(iatom.type()=="Xe") {return 2.16;}
			if(iatom.type()=="Pt") {return 1.75;}
			if(iatom.type()=="Au") {return 1.66;}
			if(iatom.type()=="Tl") {return 1.96;}
			if(iatom.type()=="Pb") {return 2.02;}
			if(iatom.type()=="U") {return 1.86;}
			else {return 1;}
		};

		int num(Atom iatom) {
			if(iatom.type()=="H") {return 1;}
			if(iatom.type()=="He") {return 2;}
			if(iatom.type()=="Li") {return 3;}
			if(iatom.type()=="Be") {return 4;}
			if(iatom.type()=="B") {return 5;}
			if(iatom.type()=="C") {return 6;}
			if(iatom.type()=="N") {return 7;}
			if(iatom.type()=="O") {return 8;}
			if(iatom.type()=="F") {return 9;}
			if(iatom.type()=="Ne") {return 10;}
			if(iatom.type()=="Na") {return 11;}
			if(iatom.type()=="Mg") {return 12;}
			if(iatom.type()=="Al") {return 13;}
			if(iatom.type()=="Si") {return 14;}
			if(iatom.type()=="P") {return 15;}
			if(iatom.type()=="S") {return 16;}
			if(iatom.type()=="Cl") {return 17;}
			if(iatom.type()=="Ar") {return 18;}
			if(iatom.type()=="K") {return 19;}
			if(iatom.type()=="Ca") {return 20;}
			else {return 1;}
		}
};

class GridSize {
	private:
		Vec m_origin;
        Vec m_delta;
        unsigned m_nx;
        unsigned m_ny;
        unsigned m_nz;

	public:
		GridSize(Vec const& origin, Vec const& delta, unsigned const nx, unsigned const ny, unsigned const nz) : 
		m_origin(origin), m_delta(delta), m_nx(nx), m_ny(ny), m_nz(nz) {}

        GridSize(Vec const& min, Vec const& max, unsigned const quality) 
        {
         	double d(stepSize(quality));

         	m_origin = min;
         	m_delta.setValue(d, d, d);

         	Vec delta(max-min);
   			delta /= d;

   			m_nx = std::ceil(delta.x) + 1;
   			m_ny = std::ceil(delta.y) + 1;
   			m_nz = std::ceil(delta.z) + 1;
         }

        Vec const& delta() const { return m_delta; }
        Vec const& origin() const { return m_origin; }

        unsigned nx() const { return m_nx; }
        unsigned ny() const { return m_ny; }
        unsigned nz() const { return m_nz; }

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
};

class GridData {
	private:

		Array3D m_data;
		Vec m_origin;
		Vec m_delta;
		
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

		Vec const& origin() const { return m_origin; }
        Vec const& delta() const {return m_delta; }

        void copy(GridData const& that)
        {
        	m_origin 	= that.m_origin;
        	m_delta 	= that.m_delta;

        	unsigned nx, ny, nz;
        	that.getNumberOfPoints(nx, ny, nz);
        	Array3D::extent_gen extents;
        	m_data.resize(extents[nx][ny][nz]);
        	m_data = that.m_data;
        }

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

		double interpolate(double const x, double const y, double const z) const {

			double weight(6.0);
			unsigned nx,ny,nz;
			getNumberOfPoints(nx,ny,nz);

			if (x==0 || y==0 || z==0) {return m_data[x][y][z];}
			if (x==(nx-1) || y==(ny-1) || z==(nz-1)) {return m_data[x][y][z];}

			double value = (1/weight)*( m_data[x-1][y][z] 
									  + m_data[x+1][y][z] 
									  + m_data[x][y-1][z]
									  + m_data[x][y+1][z]
									  + m_data[x][y][z-1]
									  + m_data[x][y][z+1] );
			return value;
		}

		bool saveToCubeFile(QString const& filePath, QStringList const& coordinates, bool const invertSign) const {
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
		double const& radius() const {return m_radius;}
};

Vec getBBMax(QVector<Atom> atomVector, Solvent solvent_instance) {
	
	double Ra = 0;
	double Rs = solvent_instance.radius();
	Vec max(0.0,0.0,0.0);

	for (int i=0; i<atomVector.size(); i++) {
		if(atomVector[i].vdw(atomVector[i]) > Ra) { Ra = atomVector[i].vdw(atomVector[i]); }
		if(atomVector[i].getX() > max.x) {max.x = atomVector[i].getX();}
		if(atomVector[i].getY() > max.y) {max.y = atomVector[i].getY();}
		if(atomVector[i].getZ() > max.z) {max.z = atomVector[i].getZ();}
	}

	Vec solvent(Rs, Rs, Rs);
	Vec vdw(Ra, Ra, Ra);
	max += solvent+vdw;

	return max;
}; 

Vec getBBMin(QVector<Atom> atomVector, Solvent solvent_instance) {

	double Ra = 0;
	double Rs = solvent_instance.radius();
	Vec min(0.0,0.0,0.0);

	for (int i=0; i<atomVector.size(); i++) {
		if(atomVector[i].vdw(atomVector[i]) > Ra) { Ra = atomVector[i].vdw(atomVector[i]); }
		if(atomVector[i].getX() < min.x) {min.x = atomVector[i].getX();}
		if(atomVector[i].getY() < min.y) {min.y = atomVector[i].getY();}
		if(atomVector[i].getZ() < min.z) {min.z = atomVector[i].getZ();}
	}

	Vec solvent(Rs, Rs, Rs);
	Vec vdw(Ra, Ra, Ra);
	min -= solvent+vdw;

	return min;
};

float DistanceBetween3DPoints(float const x1, float const x2, 
	float const y1, float const y2, float const z1, float const z2) {

	float distance = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );
	return distance;
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

GridData doVDWTest(GridSize size, QVector<Atom> atoms, GridData data) {

	for (unsigned i=0; i<size.nx(); i++) {
		for (unsigned j=0; j<size.ny(); j++) {
			for (unsigned k=0; k<size.nz(); k++) {
				
				float x1 = size.origin()[0] + i*size.delta()[0];
				float y1 = size.origin()[1] + j*size.delta()[1];
				float z1 = size.origin()[2] + k*size.delta()[2];

				for (int ia=0; ia<atoms.size(); ia++) {

					float x2 = atoms[ia].getX();
					float y2 = atoms[ia].getY();
					float z2 = atoms[ia].getZ();

					if (DistanceBetween3DPoints(x1,x2,y1,y2,z1,z2) <= atoms[ia].vdw(atoms[ia])) {
						data.operator()(i,j,k) = 1;
						break;
					}
				}
			}
		}
	}

	return data;
}

GridData doSESTest(GridSize size, QVector<Atom> atoms, GridData data, Solvent solvent) {

	double Rs = solvent.radius();
	int NPWS = Rs/size.delta()[0];
	QVector<GridPoint> admissible;

	for (unsigned i=0; i<size.nx(); i++) {
		for (unsigned j=0; j<size.ny(); j++) {
			for (unsigned k=0; k<size.nz(); k++) {

				float x1 = size.origin()[0] + i*size.delta()[0];
				float y1 = size.origin()[1] + j*size.delta()[1];
				float z1 = size.origin()[2] + k*size.delta()[2];

				// Pretend point is admissible, then invalidate this if it is not true.
				data.operator()(i,j,k) = 1;

				for (int iatom=0; iatom<atoms.size(); iatom++) {

					float x2 = atoms[iatom].getX();
					float y2 = atoms[iatom].getY();
					float z2 = atoms[iatom].getZ();

					if (DistanceBetween3DPoints(x1,x2,y1,y2,z1,z2) < (atoms[iatom].vdw(atoms[iatom]) + Rs)) {
						data.operator()(i,j,k) = 0;
						break;
					}
				}

				if (data.operator()(i,j,k) == 1) {
					GridPoint admissiblePoint(i,j,k);
					admissible.append(admissiblePoint);
				}
			}
		}
	}

	for (int indx=0; indx<admissible.size(); indx++) {

		unsigned i = get<0>(admissible[indx]);
		unsigned j = get<1>(admissible[indx]);
		unsigned k = get<2>(admissible[indx]);

		float x1 = size.origin()[0] + i*size.delta()[0];
		float y1 = size.origin()[1] + j*size.delta()[1];
		float z1 = size.origin()[2] + k*size.delta()[2];

		for (int d1 = -NPWS; d1 < NPWS+1; d1++) {
			for (int d2 = -NPWS; d2 < NPWS+1; d2++) {
				for (int d3 = -NPWS; d3 < NPWS+1; d3++) {

					unsigned int ind1 = i + d1;
					unsigned int ind2 = j + d2;
					unsigned int ind3 = k + d3;

					// Trying to index outside the grid itself (less than minimum).
					if (std::signbit(ind1) || std::signbit(ind2) || std::signbit(ind3)) {
						continue;
					}

					// Trying to index outside the grid itself (greater than maximum) - Remember zero indexing.
					if (ind1 >= size.nx() || ind2 >= size.ny() || ind3 >= size.nz()) {
						continue;
					}

					float x2 = size.origin()[0] + (ind1)*size.delta()[0];
					float y2 = size.origin()[1] + (ind2)*size.delta()[1];
					float z2 = size.origin()[2] + (ind3)*size.delta()[2];

					if ( DistanceBetween3DPoints(x1,x2,y1,y2,z1,z2) <= Rs) {
						data.operator()(ind1, ind2, ind3) = 1;
					}
				}
			}
		}
	}

	return data;
}

GridData doGridSmoothing(GridSize size, GridData data, bool gridsmoothing) {

	if (gridsmoothing) {
		// Create a copy of data to be modfied
		GridData datacopy(size);
		datacopy.copy(data);

		for (unsigned i=0; i<size.nx(); i++) {
			for (unsigned j=0; j<size.ny(); j++) {
				for (unsigned k=0; k<size.nz(); k++) {
					/// Modify the grid copy only. Use the original grid as an unmutable reference point.
					datacopy.operator()(i,j,k) = data.interpolate(i,j,k);
				}
			}
		}
		return datacopy;
	}
	else {
		return data;
	}
}

QStringList getCoordinates(QVector<Atom> atomVector) {

	QStringList coordinates;
	for (int iatom = 0; iatom < atomVector.size(); iatom++) {
		coordinates << QString("%1 %2 %3 %4 %5").arg(atomVector[iatom].num(atomVector[iatom]), 5)
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

	int Quality = 4;
	bool gridsmoothing = true;

	QVector<Atom> atomVector = parseXYZ("ring2.xyz");
	QString filePath = "ring2-SEStest-300817-smooth.cube";
	Solvent DefaultSolvent("Water", 1.4);

	Vec BBMin = getBBMin(atomVector, DefaultSolvent);
	Vec BBMax = getBBMax(atomVector, DefaultSolvent);

	GridSize myGridSize(BBMin, BBMax, Quality);
	GridData myGridData(myGridSize);	

	if (DefaultSolvent.radius()/myGridSize.stepSize(Quality) < 2) {
		cout << endl << "WARNING: You may need to increase the density of your grid!" << endl;
	}

	GridData test = doSESTest(myGridSize, atomVector, myGridData, DefaultSolvent);
	GridData smoothTest = doGridSmoothing(myGridSize, test, gridsmoothing);

	int testTimeNS = myTimer.elapsed(); // time taken in ns
	cout << " " << endl;
	cout << "SESTest completed in " << testTimeNS << " ns @ Quality = " << Quality << "." << endl;
	cout << " " << endl;

	QStringList coordinates = getCoordinates(atomVector);

	smoothTest.saveToCubeFile(filePath, coordinates, false);

	return 0;
}




