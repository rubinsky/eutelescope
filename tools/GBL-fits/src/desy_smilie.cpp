
// Igor Rubinskiy, Sep 2012 - Nov 2013
// GBL interpolation study: | | | [] | | | 30 mm spacing
// CERN: 120 GeV
// DESY: 1-2-3-4-5 GeV 
// unit: mm

#include <time.h>
#include "dut30.h" // empty
#include "TRandom3.h"
#include "GblTrajectory.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>



using namespace std;
//----------------------------------------------------------------------------

//  
    int addMaterial = 0;

    std::vector<Double_t> pbeam;
    unsigned int iLabel = 0;
    double step = 0;
    double s = 0;

    int counter=0;
    int dbg =  0;

//  int dz = 20;   // [mm] z spacing

    double resx =  3.0e-03;//3.0E-6; // [mm] intrinsic resolution cols, 00 deg
    double resy =  3.0e-03;//3.0E-6; // [mm] intrinsic resolution rows, 00 deg
 
    double tetSi     = 0.;
    double tetDUT    = 0.;
    double tetAIR    = 0.;
    double tetKapton = 0.;
    double tet2Plane = 0.;

    double X0Kapton  = 0 ; // Kapton (Polyimide Film)
    double X0Si      = 0 ; // Sensor + ROC + PCB
    double X0DUT     = 0 ;
    double X0AIR     = 0 ;
    double X02Plane  = 0 ;

    TMatrixD jacPointToPoint(5, 5);

    GblPoint *point = 0;
    TMatrixD proL2m(2,2);

    int nplanes = 7;
    int kplane  = 0;

    std::vector<double_t> iz;
 
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

//  jac[1][3] = ds;//x = xp * ds
//  jac[2][4] = ds;//y = yp * ds

  return jac;
}

double GetScatter(double p, double XoverX0  )
{
    return  0.0136*pow(p, -1.000) * pow( XoverX0, 0.500 ) * ( 1 + 0.038*1.000*log( XoverX0) ) * 1.000 ;
//    return  0.0136*pow(p, -1.000) * pow( XoverX0, 0.5555 ) ;
}

double GetScatter(double p, double thickness , double material  )
{
//      _planeScat[ipl]= 
// pow(GetBeamValue(),GetParameter(1)) * 
// pow( GetPlaneThicknessValue( ipl ) /_planeX0[ipl], GetParameter(2) ) 
//      * (1.+0.038*GetParameter(3)*std::log( GetPlaneThicknessValue( ipl ) /_planeX0[ipl])) *  GetParameter(4) ;
// 17    0.066    0.000
// 18   -1.027    0.001
// 19    0.560    0.002
// 20    3.680    0.001
// 21    0.937    0.012

   double XoverX0   = ( thickness ) / material;

    return  0.0136*pow(p, -1.000) * pow( XoverX0, 0.500 ) * ( 1 + 0.038*1.000*log( XoverX0) ) * 1.000 ;
//   return  0.0136*pow(p, -1.027) * pow( XoverX0, 0.560 ) * ( 1 + 0.038*3.680*log( XoverX0) ) * 0.937 ;
}

//----------------------------------------------------------------------------
int addKapton(GblTrajectory &traj, TVectorD &meas, TVectorD &measPrec, TVectorD &scat, TVectorD &wscatSi, TVectorD &wscatDUT, std::vector<double> &sPoint, bool dut)
{
return 0;
}

//----------------------------------------------------------------------------
int BookPlane(GblTrajectory &traj, TVectorD &meas, TVectorD &measPrec, TVectorD &scat, TVectorD &wscatDUT, std::vector<double> &sPoint, bool dut)
{

    //----------------------------------------------------------------------------
    // plane iLabel:
    s += step;//arc length

    TMatrixD jacPointToPoint(5, 5);
    jacPointToPoint.UnitMatrix();

    jacPointToPoint = Jac5( step );
    point = new GblPoint(jacPointToPoint);

    proL2m.UnitMatrix();

    if(dut == 0) 
    {
      point->addMeasurement( proL2m, meas, measPrec );
      point->addScatterer( scat, wscatDUT );
    }
    else
    {
      point->addScatterer( scat, wscatDUT );
    }

    // point->printPoint();
    // add point to trajectory:
    iLabel = traj.addPoint(*point);

    sPoint.push_back(s);

    if( dbg > 1 )
    {
      cout << "  " << iLabel;
      cout << " at " << s;
      cout << " dut " << dut;
      cout << "  hasMeas " << point->hasMeasurement() << " res " << resy*1E3;
      cout << ", hasScat " << point->hasScatterer() << " rms " << tetSi*1E3;
      cout << endl;
    }
    delete point;

    return iLabel;
}


int main(int argc, const char* argv[]) {

   if(argc > 0)
   {
     addMaterial = atoi(argv[1]); 
   }
   else
   {
     addMaterial = 1;
   }
  // user input:

  // double resx =  1.0E-4; // [mm] intrinsic resolution rows, 0 deg
  // double resy =  1.0E-4; // [mm] intrinsic resolution cols, 0 deg

  // X0 Si = 21.82/2.33 = 9.365 cm

  //----------------------------------------------------------------------------

  clock_t startTime = clock();

  // measurement = residual
  TVectorD meas(2);
  meas.Zero();//ideal

  TVectorD measPrec(2); // precision = 1/resolution^2

  // scatter:
  TVectorD scat(2);
  scat.Zero();//mean is zero

  TVectorD wscatSi(2);
  TVectorD wscatDUT(2);
  TVectorD wscatAIR(2);
  TVectorD wscatKAPTON(2);
  TVectorD wscat2Plane(2);



// loop over momentum:

//  double dzplus=10.;
/*
  pbeam.push_back(0.1);
  pbeam.push_back(1.);
  pbeam.push_back(2.);
  pbeam.push_back(3.);
  pbeam.push_back(4.);
  pbeam.push_back(5.);
  pbeam.push_back(6.);
  pbeam.push_back(15.);
  pbeam.push_back(120.);
*/
//  pbeam.push_back(1e3);

  for(double p=1000.0;p>0.1;p=p/2.)
  {
    pbeam.push_back(p);
  }




  do{
 
  measPrec[0] = 1.0 / resx / resx;
  measPrec[1] = 1.0 / resy / resy;


  std::map<int, double> pp;
  std::map<int, map<int,double> > rx;
  std::map<int, map<int,double> > ry;

  std::map<int, map<int,double> > thx;
  std::map<int, map<int,double> > thy;
  double p = 1.0;
  int ip = 0;
 
  for(double ddutth = 50.; ddutth<=50000.; ddutth*=10  )
  {

  for(double ddut = 20.; ddut<=150.; ddut+=130 )
  {

  for(double dz = 20.; dz<=150.; dz+=130  )
  {

  for(double dz2 = 20.; dz2<=150.; dz2+=130  )
  {


  X0Kapton = ( 25e-03   ) / 286  ; // Kapton (Polyimide Film)
  X0Si     = ( 50e-03   ) / 93.6 ; // Sensor + ROC + PCB
  X0DUT    = ( ddutth*1e-03   ) / 93.6 ;
  X0AIR    = ( dz /* mm */ ) / 300000 /* 300 m is rad length of the air */ ;

  iz.clear();
  //printf("iz : %5d \n", iz.size());
  double z_plane=0.;
// 0   1   2   3   4   5   6
//  dz  dz2 dut dut dz2 dz
  for(int k = 0; k < nplanes; k++)
  {
    iz.push_back(z_plane);
    if( k<1  || k >4  ) z_plane += dz;
    else if( k<2  || k >3  ) z_plane += dz2;
    if( k==2 || k==3  ) z_plane += ddut;
   // printf("dz %5.2f dz2 %5.2f ddut %5.2f k=%2d iz:%5.2f\n",dz,dz2,ddut,k,iz[k]);
  }

  for(int i=0; i<pbeam.size(); i++)
  {
    p = pbeam[i];

    if( dbg > 1)
    {
      cout << endl;
      cout << "p " << p << " GeV" << endl;
    }

    pp[i] = p;

//    tetSi  = 0.0136 * sqrt( X0Si) / p * ( 1 + 0.038*log( X0Si) );
//    tetDUT = 0.0136 * sqrt(X0DUT) / p * ( 1 + 0.038*log(X0DUT) );
//    tetAIR = 0.0136 * sqrt(X0AIR) / p * ( 1 + 0.038*log(X0AIR) );

//    tetSi     = GetScatter( p, 0.050, 93.6);
//    tetDUT    = GetScatter( p, 0.050, 93.6);
//    tetAIR    = GetScatter( p, dz , 300000. );
//    tetKapton = GetScatter( p, 0.050 , 286. );

    tetSi     = GetScatter( p, X0Si );
    tetDUT    = GetScatter( p, X0DUT );
    tetAIR    = GetScatter( p, X0AIR );
    tetKapton = GetScatter( p, X0Kapton );
    tet2Plane = GetScatter( p, (X0AIR+X0Kapton*2)/2.),// devide by 2 to have 2 thin scatterers

    wscatSi[0]  = 1.0 / ( tetSi * tetSi );//weight
    wscatSi[1]  = 1.0 / ( tetSi * tetSi );

    wscatDUT[0] = 1.0 / ( tetDUT * tetDUT );//weight
    wscatDUT[1] = 1.0 / ( tetDUT * tetDUT );

    wscatAIR[0] = 1.0 / ( tetAIR * tetAIR );//weight
    wscatAIR[1] = 1.0 / ( tetAIR * tetAIR );

    wscatKAPTON[0] = 1.0 / ( tetKapton * tetKapton );//weight
    wscatKAPTON[1] = 1.0 / ( tetKapton * tetKapton );

    wscat2Plane[0] = 1.0 / ( tet2Plane * tet2Plane );//weight
    wscat2Plane[1] = 1.0 / ( tet2Plane * tet2Plane );


    std::vector<double> sPoint;

    // create trajectory:
    printf("dz: %3.0f dz2: %3.0f ddut: %3.0f  p: %07.2f  %-4.2f um :: ", dz,dz2,ddut, p, resx*1e3);

    for(int k = 0; k < nplanes; k++)
    {
      kplane = k;
      GblTrajectory traj( false ); // curvature = false

    
      std::vector<int>  itel(nplanes); 
      std::vector<bool> idut(nplanes);
      for(int i=0;i<nplanes;i++)
      {
        itel[i] = 0; idut[i] = 0;     
        if(i==kplane) { idut[i] = 1; }
      }

//      step = dz;
 
//        double dev = 3.0;
//        if(dev>2)
//        {
//          step = dz /dev;
//          s=-step;
//        }
 
//      printf("addMaterial: %d \n", addMaterial);

      for(int i=0;i<nplanes;i++)
      {
          int dummy =0;
          if(i==0)
          {
            step = 0;
            itel[i] = BookPlane( traj, meas, measPrec, scat, wscatSi, sPoint, idut[i]);
          }
          if(i<1) continue;

          if(addMaterial == 1)
          {
          // continue if i > 0
          step = iz[i]-iz[i-1];
          //printf("z1 %5f z2 %5f step %5f \n",iz[i],iz[i-1],step);
          double step1= step/2. - step/sqrt(12); // first scatterer
          double step2= 2 * step/sqrt(12);       // second scatterer
          double step3= step1;                   // sensor

          step = step1;
          dummy = BookPlane( traj, meas, measPrec, scat, wscat2Plane, sPoint, 1);

          step = step2;
          dummy = BookPlane( traj, meas, measPrec, scat, wscat2Plane, sPoint, 1);

          step = step3;
          if(i==3)
          {
            itel[i] = BookPlane( traj, meas, measPrec, scat, wscatDUT,sPoint, idut[i]);
          }
          else if(i=nplanes-1)
          { 
             step = 1000.;
             itel[i] = BookPlane( traj, meas, measPrec, scat, wscatSi, sPoint, idut[i]);
          }
          else if(i!=3)
          { 
              itel[i] = BookPlane( traj, meas, measPrec, scat, wscatSi, sPoint, idut[i]);
          }


/*
//add Kapton
          step = ( dz-15.)/2.;
          int dummy = BookPlane( traj, meas, measPrec, scat, wscatKAPTON, sPoint, 1);

//add sensors
          step =13.0;
          itel[i] = BookPlane( traj, meas, measPrec, scat, wscatSi, sPoint, idut[i]);

//add Kapton
          step = 2.0;
          dummy = BookPlane( traj, meas, measPrec, scat, wscatKAPTON, sPoint, 1);

//add air
          step = (dz-15.)/2.;       
          dummy = BookPlane( traj, meas, measPrec, scat,  wscatAIR, sPoint, 1 );
*/
        }else{
//add sensors
          step = iz[i]-iz[i-1];
          if(i!=3) itel[i] = BookPlane( traj, meas, measPrec, scat, wscatSi, sPoint, idut[i]);
          if(i==3) itel[i] = BookPlane( traj, meas, measPrec, scat, wscatDUT,sPoint, idut[i]);
        }
      }

// addKapton( traj, meas, measPrec, scat, wscatSi, wscatDUT, sPoint, idut[i]);
// fit trajectory:

      double Chi2;
      int Ndf;
      double lostWeight;

      traj.fit( Chi2, Ndf, lostWeight );

//    cout << " Fit: " << Chi2 << ", " << Ndf << ", " << lostWeight << endl;

      TVectorD aCorrection(5);
      TMatrixDSym aCovariance(5);

// at DUT:

// results for track parameters :
      int ipos = itel[kplane];
      traj.getResults( ipos, aCorrection, aCovariance );
// Covariance matric elements: 1,1 - x' - the error on the slope (dx/dz) based on the errors, MS, geometry from other planes

// getMeasResults for points use in the fit
//      unsigned int numData;  // 2 elements X and Y      
//      traj.getMeasResults(ipos, numData,   aResiduals, aMeasErrors,  	aResErrors, aDownWeights);
// aResiduals  = 0.0 (in this case only!)
// aMeasErrors = 4.0 um (the measurement error = intrinisic Mimosa26 resolution)
// aResErrors  = biased residual RMS

// for scattering planes:
//      traj.getScatResults(ipos, numData,   aResiduals, aMeasErrors,  	aResErrors, aDownWeights);
// aResiduals  - 0 in this case
// aMeasErrors - the error on the MS (the theta I have put in with getScatter method f(X of the material) )
// aResErrors  - error of the kink

//    cout << endl;
//    cout << "point " << ipos;
//    cout << " at " << sPoint[ipos-1] << " mm:";
//    cout << "  sigma(x) = " << sqrt(aCovariance(3,3))*1E3 << " um";
//    cout << ", sigma(y) = " << sqrt(aCovariance(4,4))*1E3 << " um";
//    cout << endl;

      rx[ip][k] = sqrt(aCovariance(3,3))*1E0;
      ry[ip][k] = sqrt(aCovariance(4,4))*1E0;

      thx[ip][k] = sqrt(aCovariance(1,1))*1E3;
      thy[ip][k] = sqrt(aCovariance(2,2))*1E3;
      printf("%-7.2f",   rx[ip][k]*1e3);
    }
 
    printf("        ");
    for(int k = 0; k < nplanes; k++)
    {
      printf("%-7.2f",   sqrt(rx[ip][k]*rx[ip][k] + resx*resx)*1e3 );
    }

/*
    printf("        ");
    for(int k = 0; k < nplanes; k++)
    {
      printf("%-4.1f",   sqrt(resx*resx - rx[ip][k]*rx[ip][k] )*1e3 );
    }
*/

    printf(" %-7.2f ", ddutth);

    printf(" ");
    for(int k = 0; k < nplanes; k++)
    {
      printf("%-7.3f", thx[ip][k]  );
    }

   printf("\n");


    ip++;


  }//p loop
  }//dz2  loop
  }//dz  loop
  }//ddut  loop

  }//ddut thickness loop



/*
  cout << endl;
  cout << "  Float_t GBLp"<< counter <<"[" << ip << "] = { ";

  for( int jp = 0; jp < ip; ++jp ) {
    cout << fixed << showpoint << setprecision(2) << pp[jp];
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
*/

/*
  cout << endl;
  cout << "  Float_t GBLrxRES"<< counter <<"[" << ip << "] = { ";
  for( int jp = 0; jp < ip; ++jp ) {
    cout << fixed << showpoint << setprecision(2) << sqrt(ry[jp]*ry[jp]+(resy*1000.)*(resy*1000.));
    if( jp < ip-1 ) cout << ", ";
  }
  cout << " };" << endl;
*/

  resx += 0.0001;
  resy += 0.0001;
  counter++;
  }//res loop
  while ( resx < 0.005);

//    dz *= 2.; // obstain
//    resx -= 0.0005;
//    resy -= 0.0005;


//  resy -= 0.0001; 
//  counter++;
//  }//dzplus loop
//  while ( resx > 0.003 ) ;


  clock_t endTime = clock();
  double diff = endTime - startTime;
  double cps = CLOCKS_PER_SEC;
  cout << endl << " Time elapsed " << diff / cps << " s" << endl;

  return 0;
}
