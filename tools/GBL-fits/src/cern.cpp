
// Daniel Pitzl, Mar 2012
// GBL interpolation study: | | | [] | | | 30 mm spacing
// CERN: 120 GeV
// unit: mm

#include <time.h>
#include "dut30.h" // empty
#include "TRandom3.h"
#include "GblTrajectory.h"
#include <cmath>
#include <iomanip>

using namespace std;

//----------------------------------------------------------------------------

TMatrixD Jac5( double ds ) {
  /*
    Jacobian for straight line track
    track = q/p, x', y', x, y
            0,   1,  2,  3, 4
  */
  TMatrixD jac(5, 5);
  jac.UnitMatrix();
  jac[3][1] = ds;//x = xp * ds
  jac[4][2] = ds;//y = yp * ds
  return jac;
}

//----------------------------------------------------------------------------

int main() {

  // user input:

  double dz =  10;// [mm] z spacing

  double resx =  4.0E-3; // [mm] intrinsic resolution cols, 00 deg
  double resy =  2.0E-3; // [mm] intrinsic resolution rows, 00 deg

  //double resx = 45.0E-3; // [mm] intrinsic resolution rows, 0 deg
  //double resy = 25.0E-3; // [mm] intrinsic resolution cols, 0 deg

  // X0 Si = 21.82/2.33 = 9.365 cm

  double X0Si = 0.05 / 94; // Sensor + ROC + PCB

  double X0DUT = 0.05;

  //----------------------------------------------------------------------------

  clock_t startTime = clock();

  // measurement = residual
  TVectorD meas(2);
  meas.Zero();//ideal

  TVectorD measPrec(2); // precision = 1/resolution^2
  measPrec[0] = 1.0 / resx / resx;
  measPrec[1] = 1.0 / resy / resy;

  // scatter:
  TVectorD scat(2);
  scat.Zero();//mean is zero

  TVectorD wscatSi(2);
  TVectorD wscatDUT(2);

  // loop over momentum:

  int counter=0;

  do{
 
  double pp[99];
  double rx[99];
  double ry[99];
  double p = 1.0;
  int ip = 0;
 
  do{

//    cout << endl;
//    cout << "p " << p << " GeV" << endl;

    pp[ip] = p;

    double tetSi = 0.0136 * sqrt(X0Si) / p * ( 1 + 0.038*log(X0Si) );

    wscatSi[0] = 1.0 / ( tetSi * tetSi );//weight
    wscatSi[1] = 1.0 / ( tetSi * tetSi );

    double tetDUT = 0.0136 * sqrt(X0DUT) / p * ( 1 + 0.038*log(X0DUT) );

    wscatDUT[0] = 1.0 / ( tetDUT * tetDUT );//weight
    wscatDUT[1] = 1.0 / ( tetDUT * tetDUT );

    std::vector<double> sPoint;

    // create trajectory:

    GblTrajectory traj( false ); // curvature = false

    double step = 0;
    double s = 0;

    //----------------------------------------------------------------------------
    // plane 0:

    TMatrixD jacPointToPoint(5, 5);
    jacPointToPoint.UnitMatrix();

    GblPoint *point = new GblPoint(jacPointToPoint);

    TMatrixD proL2m(2,2);
    proL2m.UnitMatrix();

    point->addMeasurement( proL2m, meas, measPrec );

    point->addScatterer( scat, wscatSi );

    // add point to trajectory:
    unsigned int iLabel = traj.addPoint(*point);

    sPoint.push_back(s);

//    cout << "  " << iLabel;
//    cout << " at " << s;
//    cout << "  hasMeas " << point->hasMeasurement() << " res " << resy*1E3;
//    cout << ", hasScat " << point->hasScatterer() << " rms " << tetSi*1E3;
//    cout << endl;

    delete point;

    //----------------------------------------------------------------------------
    // plane 1:

    step = dz;
    jacPointToPoint = Jac5( step );
    point = new GblPoint(jacPointToPoint);
    point->addMeasurement( proL2m, meas, measPrec );
    point->addScatterer( scat, wscatSi );
    iLabel = traj.addPoint(*point);

    s += step;//arc length
    sPoint.push_back(s);

//    cout << "  " << iLabel;
//    cout << " at " << s;
//    cout << "  hasMeas " << point->hasMeasurement() << " res " << resy*1E3;
//    cout << ", hasScat " << point->hasScatterer() << " rms " << tetSi*1E3;
//    cout << endl;

    delete point;

    //----------------------------------------------------------------------------
    // plane 2:

    step = 10.0;
    jacPointToPoint = Jac5( step );
    point = new GblPoint(jacPointToPoint);
    point->addMeasurement( proL2m, meas, measPrec );
    point->addScatterer( scat, wscatSi );

    iLabel = traj.addPoint(*point);

    s += step;//arc length
    sPoint.push_back(s);

//    cout << "  " << iLabel;
//    cout << " at " << s;
//    cout << "  hasMeas " << point->hasMeasurement() << " res " << resy*1E3;
//    cout << ", hasScat " << point->hasScatterer() << " rms " << tetSi*1E3;
//    cout << endl;

    delete point;

    // DUT:

    step = 10.0;
    jacPointToPoint = Jac5( step );
    point = new GblPoint(jacPointToPoint);
    point->addScatterer( scat, wscatDUT );// DUT
    iLabel = traj.addPoint(*point);
    int idut = iLabel;
    s += step;//arc length
    sPoint.push_back(s);

//    cout << "  " << iLabel;
//    cout << " at " << s;
//    cout << ", hasScat " << point->hasScatterer() << " rms " << tetDUT*1E3;
//    cout << endl;

    delete point;

    //----------------------------------------------------------------------------
    // plane 3:

    step = 10.0;
    jacPointToPoint = Jac5( step );
    point = new GblPoint(jacPointToPoint);
    point->addMeasurement( proL2m, meas, measPrec );
    point->addScatterer( scat, wscatSi );
    iLabel = traj.addPoint(*point);
    s += step;//arc length
    sPoint.push_back(s);

//    cout << "  " << iLabel;
//    cout << " at " << s;
//    cout << "  hasMeas " << point->hasMeasurement();
//    cout << "  hasScat " << point->hasScatterer() << " rms " << tetSi*1E3;
//    cout << endl;

    delete point;

    //----------------------------------------------------------------------------
    // plane 4:

    step = 10.0;
    jacPointToPoint = Jac5( step );
    point = new GblPoint(jacPointToPoint);
    point->addMeasurement( proL2m, meas, measPrec );
    point->addScatterer( scat, wscatSi );
    iLabel = traj.addPoint(*point);
    s += step;//arc length
    sPoint.push_back(s);

//    cout << "  " << iLabel;
//    cout << " at " << s;
//    cout << "  hasMeas " << point->hasMeasurement();
//    cout << "  hasScat " << point->hasScatterer() << " rms " << tetSi*1E3;
//    cout << endl;

    delete point;

    //----------------------------------------------------------------------------
    // plane 5:

    step = dz;
    jacPointToPoint = Jac5( step );
    point = new GblPoint(jacPointToPoint);
    point->addMeasurement( proL2m, meas, measPrec );
    point->addScatterer( scat, wscatSi );
    iLabel = traj.addPoint(*point);
    s += step;//arc length
    sPoint.push_back(s);

//    cout << "  " << iLabel;
//    cout << " at " << s;
//    cout << "  hasMeas " << point->hasMeasurement();
//    cout << "  hasScat " << point->hasScatterer() << " rms " << tetSi*1E3;
//    cout << endl;

    delete point;

    //----------------------------------------------------------------------------
    // fit trajectory:

    double Chi2;
    int Ndf;
    double lostWeight;

    traj.fit( Chi2, Ndf, lostWeight );

    //cout << " Fit: " << Chi2 << ", " << Ndf << ", " << lostWeight << endl;

    TVectorD aCorrection(5);
    TMatrixDSym aCovariance(5);

    // at DUT:

    int ipos = idut;
    traj.getResults( ipos, aCorrection, aCovariance );

//    cout << endl;
//    cout << "point " << ipos;
//    cout << " at " << sPoint[ipos-1] << " mm:";
//    cout << "  sigma(x) = " << sqrt(aCovariance(3,3))*1E3 << " um";
//    cout << ", sigma(y) = " << sqrt(aCovariance(4,4))*1E3 << " um";
//    cout << endl;

    rx[ip] = sqrt(aCovariance(3,3))*1E3;
    ry[ip] = sqrt(aCovariance(4,4))*1E3;

    ip++;

    p *= 1.1;

  }//p loop
  while ( p < 1000 );


  cout << endl;
  cout << "  Float_t GBLp"<< counter <<"[" << ip << "] = { ";

  for( int jp = 0; jp < ip; ++jp ) {
    cout << fixed << showpoint << setprecision(3) << pp[jp];
    if( jp < ip-1 ) cout << ", ";
  }
  cout << " };" << endl;

  cout << endl;
  cout << "  Float_t GBLrx"<< counter <<"[" << ip << "] = { ";
  for( int jp = 0; jp < ip; ++jp ) {
    cout << fixed << showpoint << setprecision(2) << rx[jp];
    if( jp < ip-1 ) cout << ", ";
  }
  cout << " };" << endl;

  cout << endl;
  cout << "  Float_t GBLry"<< counter <<"[" << ip << "] = { ";
  for( int jp = 0; jp < ip; ++jp ) {
    cout << fixed << showpoint << setprecision(2) << ry[jp];
    if( jp < ip-1 ) cout << ", ";
  }
  cout << " };" << endl;


  dz *= 2.; counter++;
  }//dz loop
  while ( dz < 330 );


  clock_t endTime = clock();
  double diff = endTime - startTime;
  double cps = CLOCKS_PER_SEC;
  cout << endl << " Time elapsed " << diff / cps << " s" << endl;

  return 0;
}
