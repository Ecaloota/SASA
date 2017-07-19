/* Parse an input .xyz file and return a .cube file containing
solvent-accessible surface area data encoded within a 
binary array, for visualisation within IQMol. */

#include <fstream>
#include "boost/multi_array.hpp"
#include <QVector>
#include <QString>
#include <QStringList>
#include <QFile>
#include <QTextStream>

//#include "QGLViewer/vec.h"

using namespace std;

typedef boost::multi_array<double, 3> Array3D;

double const BohrRadius          = 5.2917721092e-11;
double const BohrToAngstrom      = BohrRadius*1.0e10;
double const AngstromToBohr      = 1.0/BohrToAngstrom;

// An atom contains only an atomType and a set of coords.
class Atom {
	private:
		string m_atomtype;
		float m_xcoord;
		float m_ycoord;
		float m_zcoord;
	public:
		
		Atom(string &atomtype, float xcoord, float ycoord, float zcoord) :
			m_atomtype(atomtype),m_xcoord(xcoord), m_ycoord(ycoord), m_zcoord(zcoord) {}

		Atom() : m_atomtype(), m_xcoord(), m_ycoord(), m_zcoord() {}

		string getAtomType() {return m_atomtype;}
		float getX() {return m_xcoord;}
		float getY() {return m_ycoord;}
		float getZ() {return m_zcoord;}

		float getvdW_radius(Atom iatom) {
			if(iatom.getAtomType()=="H") {return 1.20;}
			if(iatom.getAtomType()=="C") {return 1.70;}
			if(iatom.getAtomType()=="O") {return 1.52;}
			else {return 1;}
		};
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

QVector<double> get_bbmax(QVector<Atom> atomVector) {
	
	float maxX = 0;
	float maxY = 0;
	float maxZ = 0;	

	for(int iatom=0; iatom < atomVector.size(); iatom++) {

		if(atomVector[iatom].getX() > maxX) {maxX = atomVector[iatom].getX();}
		if(atomVector[iatom].getY() > maxY) {maxY = atomVector[iatom].getY();}
		if(atomVector[iatom].getZ() > maxZ) {maxZ = atomVector[iatom].getZ();}
	}
	
	maxX += 2; // fudge factor
	maxY += 2; // fudge factor
	maxZ += 2; // fudge factor

	QVector<double> bbmax = {maxX, maxY, maxZ};

	return bbmax;
}; 

QVector<double> get_bbmin(QVector<Atom> atomVector) {

	float minX = 0;
	float minY = 0;
	float minZ = 0;

	for(int iatom=0; iatom < atomVector.size(); iatom++) {
		if(atomVector[iatom].getX() < minX){minX = atomVector[iatom].getX();}
		if(atomVector[iatom].getY() < minY){minY = atomVector[iatom].getY();}
		if(atomVector[iatom].getZ() < minZ){minZ = atomVector[iatom].getZ();}
	}

	minX -= 2; // fudge factor
	minY -= 2; // fudge factor
	minZ -= 2; // fudge factor

	QVector<double> bbmin = {minX, minY, minZ};

	return bbmin;
};

QVector<Atom> parseXYZ( const string& infile ) {
	ifstream input_file(infile.c_str(), std::ifstream::in);
	if(!input_file.good()) {throw "The file does not exist! ";}

	unsigned natoms;
	float xcoord, ycoord, zcoord;
	string atomtype; 
	QVector<Atom> atomVector;

	input_file >> natoms;
	for(unsigned iatom = 0; iatom < natoms; iatom++) {
		input_file >> atomtype >> xcoord >> ycoord >> zcoord;
		atomVector.append(Atom(atomtype, xcoord, ycoord, zcoord));
	};
	
	return atomVector;
}

int main( ) 
{
	int quality = 1;
	QVector<Atom> atomVector = parseXYZ("methanol.xyz");
	QVector<double> bbmin = get_bbmin(atomVector);
	QVector<double> bbmax = get_bbmax(atomVector);

	GridSize mygridsize(bbmin, bbmax, quality);
	GridData mygriddata(mygridsize);

	for (unsigned i=0; i<mygridsize.nx(); i++) {
		for (unsigned j=0; j<mygridsize.ny(); j++) {
			for (unsigned k=0; k<mygridsize.nz(); k++) {
				double gridpoint[3] = {(mygridsize.origin()[0] + i*mygridsize.delta()[0]), (mygridsize.origin()[1] + j*mygridsize.delta()[1]), (mygridsize.origin()[2] + k*mygridsize.delta()[2])};
				for (int iatom=0; iatom<atomVector.size(); iatom++) {
					if ( sqrt( (pow((gridpoint[0]-atomVector[iatom].getX()),2)) + (pow((gridpoint[1]-atomVector[iatom].getY()),2)) + (pow((gridpoint[2]-atomVector[iatom].getZ()),2))) < atomVector[iatom].getvdW_radius(atomVector[iatom])) {
						mygriddata.operator()(i,j,k) = 1;
						break;
					}
				}
			}
		}
	}

	QStringList coordinates;
	coordinates << "1	0.00000		0.00000		1.03027		1.04441";
	coordinates << "6   0.00000		0.00000     -0.01241    0.71501";
	coordinates << "1	0.00000		0.89335 	-0.51771    1.10231";
	coordinates << "1 	0.00000		-0.89335    -0.51771    1.10231";
	coordinates << "8   0.00000 	0.00000		0.06173		-0.67782";
	coordinates << "1   0.00000     0.00000     -0.82687    -1.01010";

	QString filePath = "methanol.cube";

	mygriddata.saveToCubeFile(filePath, coordinates, false);

	return 0;
}




