#include "TString.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMinuit.h"
#include "TString.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMath.h"

#include <sstream>

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()


// system includes <>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <cstdlib>
#include <stdarg.h>
#include <limits>
#include <utility>
 #include <iomanip> 
 TCanvas *c1;
 TCanvas *c2;
 TCanvas *c3;
 TCanvas *c4;
 TH1D    *h1;
 TH1D    *h2;

 TH1D    *h3;
 TH1D    *h4;

 TLegend *leg;
 TLegend *lsmil;

 TLine *iline02; 
 TLine *iline03; 
 TLine *iline05; 
 TLine *iline15;
 TLine *iline120;

 TLine *kline02; 
 TLine *kline03; 
 TLine *kline05; 
 TLine *kline15;
 TLine *kline120;
 TGraphErrors *g;
 TGraphErrors *h;


 Int_t ipl,jc1,jc2,jc3;
// [dz01][dz12][dz23=ddut] [ibeam] [ires] [ddutth] [ipl] -> dtel

 typedef std::map<Int_t, Double_t>       t_dkapton_dtel;
 typedef std::map<Int_t, t_dkapton_dtel>   t_ipl_dkapton;
 typedef std::map<Int_t, t_ipl_dkapton>    t_dsith_ipl;
 typedef std::map<Int_t, t_dsith_ipl>    t_ddutth_dsith;
 typedef std::map<Int_t, t_ddutth_dsith> t_ires_ddutth;
 typedef std::map<Int_t, t_ires_ddutth>  t_ibeam_ires;
 typedef std::map<Int_t, t_ibeam_ires>   t_ddut_ibeam;
 typedef std::map<Int_t, t_ddut_ibeam>   t_dz2_ddut;
 typedef std::map<Int_t, t_dz2_ddut>     t_dz_dz2;


 typedef t_dkapton_dtel::iterator          t_dkapton_dtel_iter;
 typedef t_ipl_dkapton::iterator           t_ipl_dkapton_iter;
 typedef t_dsith_ipl::iterator           t_dsith_ipl_iter;
 typedef t_ddutth_dsith::iterator        t_ddutth_dsith_iter;
 typedef t_ires_ddutth::iterator         t_ires_ddutth_iter;
 typedef t_ibeam_ires::iterator          t_ibeam_ires_iter;
 typedef t_ddut_ibeam::iterator          t_ddut_ibeam_iter;
 typedef t_dz2_ddut::iterator            t_dz2_ddut_iter;
 typedef t_dz_dz2::iterator              t_dz_dz2_iter;

 t_dz_dz2 pointresolution;
 t_dz_dz2 scattering;

 std::vector<Double_t> vbeam;
 std::vector<Double_t> rbeam;


//   dz            dz2            ddut           res           thickness
   std::map<int,  std::map<int, std::map<int, std::map<int,  std::map<int,  std::map<int,  std::map<int,   std::map<int, TGraphErrors*> > > > > > > > TGres;
   std::map<int,  std::map<int, std::map<int, std::map<int,  std::map<int,  std::map<int,  std::map<int,   std::map<int, TGraphErrors*> > > > > > > > TGscat;

   std::map<int,  std::map<int, std::map<int, std::map<int,  std::map<int,  std::map<int,  std::map<int,   std::map<int, TGraphErrors*> > > > > > > > TSres;
   std::map<int,  std::map<int, std::map<int, std::map<int,  std::map<int,  std::map<int,  std::map<int,   std::map<int, TGraphErrors*> > > > > > > > TSscat;


//  dz            beam          threshold
 typedef std::map<Int_t, Double_t>     t_plane_meas;
 typedef std::map<Int_t, t_plane_meas> t_thr_plane;
 typedef std::map<Int_t, t_thr_plane>  t_beam_thr;
 typedef std::map<Int_t, t_beam_thr>   t_dz_beam;
 t_dz_beam                             measurement;
 t_dz_beam                             ermeasurement;


 typedef t_plane_meas::iterator t_plane_meas_iter;
 typedef t_thr_plane::iterator  t_thr_plane_iter;
 typedef t_beam_thr::iterator   t_beam_thr_iter;
 typedef t_dz_beam::iterator    t_dz_beam_iter;

 std::map<int, std::map<int,  std::map<int,  TGraphErrors*> > >   TGmeas;
 std::map<int, std::map<int,  std::map<int,  TGraphErrors*> > >   TSmeas;
 
 std::map<int, Double_t>  vthr;

#define nplanes 6

std::vector<std::string> parse(std::string s, std::string delimiter)
{

size_t pos = 0;
std::vector<std::string> token;
while ((pos = s.find(delimiter)) != std::string::npos) {
    token.push_back( s.substr(0, pos) );
    s.erase(0, pos + delimiter.length());
}
    token.push_back( s.substr(0, pos) );

return  token;
}

void read_DATURA()
{


//  1   2       3   4     5 6   7       8     9     10     11    12    13    14     7+8      
//   20 beam:  2.00 thr:  4 ::  0[07.24 00.05 07.32 00.06  2.016 0.008 2.008 0.008] 1[04.22 00.03 04.35 00.03  2.048 0.008 2.055 0.008] 2[04.48 00.03 04.51 00.03  2.029 0.008 2.026 0.008] 3[04.52 00.03 04.52 00.03  2.016 0.008 2.025 0.008] 4[04.44 00.03 04.48 00.03  1.999 0.008 2.013 0.008] 5[07.23 00.06 07.26 00.05  2.002 0.008 2.005 0.008]  0 00.00 00.00 00.00 00.00     1 00.00 00.00 00.00 00.00     2 00.00 00.00 00.00 00.00     3 00.00 00.00 00.00 00.00     4 00.00 00.00 00.00 00.00     5 00.00 00.00 00.00 00.00   

  Int_t idz,iplane[nplanes],ithr;
  Double_t dbeam,dthr,unres[nplanes],unerr[nplanes];
  std::string splane[nplanes];

  std::string line;
  std::string s;
  std::ifstream in("data_DATURA/A.txt",std::ios_base::in);
  
  while( std::getline(in, line) )
{
std::istringstream iss(line);

if(
iss
>> idz 
>> s >> dbeam
>> s >> ithr 
>> s 
>> splane[0] >> unerr[0] >>  s >> s >>  s >> s >> s >> s 
>> splane[1] >> unerr[1] >>  s >> s >>  s >> s >> s >> s 
>> splane[2] >> unerr[2] >>  s >> s >>  s >> s >> s >> s 
>> splane[3] >> unerr[3] >>  s >> s >>  s >> s >> s >> s 
>> splane[4] >> unerr[4] >>  s >> s >>  s >> s >> s >> s 
>> splane[5] >> unerr[5] >>  s >> s >>  s >> s >> s >> s 
){
std::vector<std::string> s;
for(int isensor = 0; isensor<nplanes; isensor++)
{
  s = parse(splane[isensor],"["); iplane[isensor]=atoi(s[0].c_str()); unres[isensor]=atof(s[1].c_str());
}

Int_t ibeam = static_cast<int> (dbeam*1.2*10);

if(vbeam.size()==0)vbeam.push_back(ibeam);
for(int i =0;i<vbeam.size();i++)
{
if(std::find(vbeam.begin(), vbeam.end(), ibeam)==vbeam.end())vbeam.push_back(ibeam);
}

printf(" %5d %5f [%5d] %5d :  %5d %-5.2f %5d %-5.2f %5d %-5.2f %5d %-5.2f %5d %-5.2f %5d %-5.2f\n",
         idz, dbeam, ibeam, ithr,
           iplane[0], unres[0], iplane[1], unres[1], iplane[2], unres[2], iplane[3], unres[3], iplane[4], unres[4], iplane[5], unres[5] );

for(int isensor = 0; isensor<nplanes; isensor++)
{
  measurement[idz][ibeam][ithr][iplane[isensor]] = unres[isensor]; ermeasurement[idz][ibeam][ithr][iplane[isensor]] = unerr[isensor];
}

}

}
 in.close();

  for( t_dz_beam_iter i0  = measurement.begin();
                      i0 != measurement.end();
                      i0++)
  { 
    for( t_beam_thr_iter i1  = i0->second.begin();
                    i1 != i0->second.end();
                    i1++)
    {
      for( t_thr_plane_iter i2  = i1->second.begin();
                      i2 != i1->second.end();
                      i2++)
      {
        for( t_plane_meas_iter i3  = i2->second.begin();
                      i3 != i2->second.end();
                      i3++)
          {
                 Int_t idz   = i0->first;
                 Int_t ibeam = i1->first;
                 Int_t ithr  = i2->first; 
                 Int_t isensor   = i3->first; 
//         std::cout  << i0->first  << " " << i1->first  << " " <<  i2->first  << " " << i3->first  << " " << i3->second << std::endl;
                  if(
                     TGmeas.find(idz) != TGmeas.end() &&                 
                     TGmeas[idz].find(ithr) != TGmeas[idz].end() &&
                     TGmeas[idz][ithr].find(isensor) != TGmeas[idz][ithr].end() 
                    )                 
                  {
//                     g = TGres[i1->first][i2->first][i4->first];
                  }
                  else
                  {  
                    TGmeas[idz][ithr][isensor] = new TGraphErrors();
                  }  
                  g = TGmeas[idz][ithr][isensor]; 
             
                  Double_t value  = i3->second;
                  Double_t value2 = value*value; 
                  Double_t vthr1 = 0.;//4.0;//vthr[i2->second]  ;
                  Double_t vthr2 = vthr1*vthr1      ;
                  Double_t pres = sqrt( value2 - vthr2 );
                  Int_t ipoint = g->GetN(); 
                  g->SetPoint( ipoint, (ibeam)/10., pres ) ;
                  g->SetPointError( ipoint, (ibeam)/10.*0.05, ermeasurement[idz][ibeam][ithr][isensor] ) ;
                  g->SetLineWidth(3);
 
                  if(
                     TSmeas.find(idz) != TSmeas.end() &&                 
                     TSmeas[idz].find(ithr) != TSmeas[idz].end() &&
                     TSmeas[idz][ithr].find(ibeam) != TSmeas[idz][ithr].end() 
                    )                 
                  {
//                     g = TGres[i1->first][i2->first][i4->first];
                  }
                  else
                  {  
                    TSmeas[idz][ithr][ibeam] = new TGraphErrors();
                  } 
                  h = TSmeas[idz][ithr][ibeam]; 
                  Int_t hpoint = h->GetN(); 
//                  std::cout<< " hpoint: " << idz << " " << ithr << " " << ibeam << " " << (isensor*idz)*1. << " " << pres << std::endl; 
                  h->SetPoint( hpoint, (isensor*idz)*1., pres ) ;
                  h->SetPointError( hpoint, 1.0, ermeasurement[idz][ibeam][ithr][isensor] ) ;
                  h->SetLineWidth(3);
//                  h->Print();
          }
      }
    }
  }

}

void read_gbl()
{
// return 0;
  //      dz            ddut           ibeam          res
  //  std::map<int, std::map<int, std::map<int, std::map< int, double > > pointresolution;
 
  std::string line;

//  TString filename("log_1_dzdz2ddut_thx.txt");
//  TString filename("log_1_dzdz2ddut_thx_6planes_kaptonAIR_stepKapton.txt");
//  TString filename("log_1_dzdz2ddut_thx_6planes_kaptonAIR_stepKapton_555.txt");
//   TString filename("log_scattering_6planes.txt");
//  TString filename("log_scattering_7planes.txt");
  TString filename("desy_smilie6_plots/resolution_6planes.txt");

  std::cout << "reading file " << filename.Data() << std::endl;
 
  // ifstream fin; 
  // fin.open(filename, ofstream::app);

  char sdummy[5];
  char sddut[4];
  char sdz[4];
  char sp[8];
  char sres[4];
  char stel[8][7];
  char sbiased[8][7];
  char sddutth[8];
 
  Int_t dz,dz2,ddut;
  Double_t pbeam,res,tel[nplanes],biased[nplanes], ddutth, dsith, dkapton;
  Double_t thx[nplanes];
  
  std::string s;
  char buf[100];
  std::ifstream in(filename.Data(),std::ios_base::in);


while( std::getline(in, line) )
{
std::istringstream iss(line);

if(iss >> s >> dz >> s >> dz2 >> s >> ddut >> s >> pbeam >> res >> s >> s 
>> tel[0] 
>> tel[1] 
>> tel[2] 
>> tel[3] 
>> tel[4] 
>> tel[5] 
//>> tel[6] 
>> biased[0] 
>> biased[1] 
>> biased[2] 
>> biased[3] 
>> biased[4] 
>> biased[5] 
//>> biased[6] 
>> ddutth
//>> thx[0]
>> thx[1]
>> thx[2]
>> thx[3]
>> thx[4]
//>> thx[5]
>> dsith
>> dkapton
){


if(dz == dz2 && dz2 == ddut )
{
printf( "%3d %3d %3d  p:%7.2f %3.1f ", dz, dz2, ddut, pbeam, res); 
printf( " tel0:%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f ", tel[0],tel[1],tel[2],tel[3],tel[4],tel[5],tel[6] );
printf( " tel0:%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f ", biased[0],biased[1],biased[2],biased[3],biased[4],biased[5],biased[6] );
//printf( " dutth: %3.1f ", ddutth);
//printf( " th1: %5.2f %5.2f %5.2f %5.2f %5.2f ", thx[1], thx[2], thx[3], thx[4], thx[5] );
//printf( " sith: %3.1f kapton:%5.1f",dsith, dkapton  );  
std::cout<<std::endl;
}
else
{
continue;
}
int ibeam  = static_cast<int>(pbeam*10);
int ires   = static_cast<int>(res*10);
int idutth = static_cast<int>(ddutth);
int isith  = static_cast<int>(dsith);
int idkapton = static_cast<int>(dkapton*10);


for(int isensor=0;isensor<nplanes;isensor++)
{
 pointresolution[dz][dz2][ddut][ibeam][ires][idutth][isith][isensor][idkapton]  = biased[isensor];//tel[3];
      scattering[dz][dz2][ddut][ibeam][ires][idutth][isith][isensor][idkapton]  = thx[isensor];
}

}else{
std::cout << "end of file" << std::endl;
}
}

// loop through 
   for( t_dz_dz2_iter i0 = pointresolution.begin();
         i0 != pointresolution.end();
         i0++)
       {
    for( t_dz2_ddut_iter i1 = i0->second.begin();
         i1 != i0->second.end();
         i1++)
       {
       for( t_ddut_ibeam_iter i2 = i1->second.begin();
            i2 != i1->second.end();
            i2++)
          {
          for( t_ibeam_ires_iter i3 = i2->second.begin();
               i3 != i2->second.end();
               i3++)
             {
              if(rbeam.size()==0)rbeam.push_back(i3->first);
              for(int i =0;i<rbeam.size();i++)
              {
               if(std::find(rbeam.begin(), rbeam.end(), i3->first )==rbeam.end())rbeam.push_back(i3->first );
              }

                for( t_ires_ddutth_iter i4 = i3->second.begin();
                  i4 != i3->second.end();
                  i4++)
                 {
// loop through ddutth
                 for( t_ddutth_dsith_iter i5 = i4->second.begin();
                     i5 != i4->second.end();
                     i5++)
                   {
                 for( t_dsith_ipl_iter i6 = i5->second.begin();
                     i6 != i5->second.end();
                     i6++)
                   {
                  for( t_ipl_dkapton_iter i7 = i6->second.begin();
                     i7 != i6->second.end();
                     i7++)
                   {
                  for( t_dkapton_dtel_iter i8 = i7->second.begin();
                     i8 != i7->second.end();
                     i8++)
                   {
/*                  std::cout  << " i0:" << i0->first  
                             << " i1 " << i1->first 
                             << " i2 " << i2->first
                             << " i3 " << i3->first 
                             << " i4 " << i4->first
                             << " i5 " << i5->first
                             << " i6 " << i6->first
                             << " i7 " << i7->first
                             << " i8 " << i8->first
                             << " i8:sec: " << i8->second << std::endl;
*/
// fill resolution .vs. momentum map
                    if(
                     TGres.find(i0->first) != TGres.end() &&                 
                     TGres[i0->first].find(i1->first) != TGres[i0->first].end() &&                 
                     TGres[i0->first][i1->first].find(i2->first) != TGres[i0->first][i1->first].end() &&
                     TGres[i0->first][i1->first][i2->first].find(i4->first) != TGres[i0->first][i1->first][i2->first].end() &&
                     TGres[i0->first][i1->first][i2->first][i4->first].find(i5->first) != TGres[i0->first][i1->first][i2->first][i4->first].end() &&
                     TGres[i0->first][i1->first][i2->first][i4->first][i5->first].find(i6->first) != TGres[i0->first][i1->first][i2->first][i4->first][i5->first].end() &&
                     TGres[i0->first][i1->first][i2->first][i4->first][i5->first][i6->first].find(i7->first) != TGres[i0->first][i1->first][i2->first][i4->first][i5->first][i6->first].end() &&
                    TGres[i0->first][i1->first][i2->first][i4->first][i5->first][i6->first][i7->first].find(i8->first) != TGres[i0->first][i1->first][i2->first][i4->first][i5->first][i6->first][i7->first].end()
                    )                 
                    {
//                     g = TGres[i1->first][i2->first][i4->first];
                    }
                    else
                    {  
                     TGres[i0->first][i1->first][i2->first][i4->first][i5->first][i6->first][i7->first][i8->first] = new TGraphErrors();
                    }  
                    g = TGres[i0->first][i1->first][i2->first][i4->first][i5->first][i6->first][i7->first][i8->first]; 
                    g->SetPoint(g->GetN(), (i3->first)/10., (i8->second)/1. );
                    g->SetLineWidth(3);
 
// fill resolution .vs. position map
                    if(
                     TSres.find(i0->first) != TSres.end() &&                 
                     TSres[i0->first].find(i1->first) != TSres[i0->first].end() &&                 
                     TSres[i0->first][i1->first].find(i2->first) != TSres[i0->first][i1->first].end() &&
                     TSres[i0->first][i1->first][i2->first].find(i4->first) != TSres[i0->first][i1->first][i2->first].end() &&
                     TSres[i0->first][i1->first][i2->first][i4->first].find(i5->first) != TSres[i0->first][i1->first][i2->first][i4->first].end() &&
                     TSres[i0->first][i1->first][i2->first][i4->first][i5->first].find(i6->first) != TSres[i0->first][i1->first][i2->first][i4->first][i5->first].end() &&
                     TSres[i0->first][i1->first][i2->first][i4->first][i5->first][i6->first].find(i3->first) != TSres[i0->first][i1->first][i2->first][i4->first][i5->first][i6->first].end() &&
                     TSres[i0->first][i1->first][i2->first][i4->first][i5->first][i6->first][i3->first].find(i8->first) != TSres[i0->first][i1->first][i2->first][i4->first][i5->first][i6->first][i3->first].end()
                    )                 
                    {
//                     g = TSres[i1->first][i2->first][i4->first];
                    }
                    else
                    {  
                     TSres[i0->first][i1->first][i2->first][i4->first][i5->first][i6->first][i3->first][i8->first]= new TGraphErrors();
                    }  
                    h = TSres[i0->first][i1->first][i2->first][i4->first][i5->first][i6->first][i3->first][i8->first]; // smilie last element is pbeam
                    h->SetPoint(h->GetN(), (i0->first * i7->first)*1., (i8->second)/1. );
                    h->SetLineWidth(3);
//h->Print();
                   }  // i8
                   } //i7
                  } //i6
//return 0;
                } //i5 
             }    //i4
          } //i3
       } //i2
    } //i1
   } //i0

// loop through 
   for( t_dz_dz2_iter i0 = scattering.begin();
         i0 != scattering.end();
         i0++)
       {
    for( t_dz2_ddut_iter i1 = i0->second.begin();
         i1 != i0->second.end();
         i1++)
       {
       for( t_ddut_ibeam_iter i2 = i1->second.begin();
            i2 != i1->second.end();
            i2++)
          {
          for( t_ibeam_ires_iter i3 = i2->second.begin();
               i3 != i2->second.end();
               i3++)
             {
             for( t_ires_ddutth_iter i4 = i3->second.begin();
                  i4 != i3->second.end();
                  i4++)
                {
// loop through ddutth
                for( t_ddutth_dsith_iter i5 = i4->second.begin();
                     i5 != i4->second.end();
                     i5++)
                   {
                 for( t_dsith_ipl_iter i6 = i5->second.begin();
                     i6 != i5->second.end();
                     i6++)
                   {
                 for( t_ipl_dkapton_iter i7 = i6->second.begin();
                     i7 != i6->second.end();
                     i7++)
                   {
                 for( t_dkapton_dtel_iter i8 = i7->second.begin();
                     i8 != i7->second.end();
                     i8++)
                   {
/*                  std::cout  << " i0:" << i0->first  
                             << " i1 " << i1->first 
                             << " i2 " << i2->first
                             << " i3 " << i3->first 
                             << " i4 " << i4->first
                             << " i5 " << i5->first
                             << " i6 " << i6->first
                             << " i7 " << i7->first
                             << " i8 " << i8->first
                             << " i8:sec: " << i8->second << std::endl;
*/
 // fill scattering.vs.momentum map
                  if(
                     TGscat.find(i0->first) != TGscat.end() &&                 
                     TGscat[i0->first].find(i1->first) != TGscat[i0->first].end() &&                 
                     TGscat[i0->first][i1->first].find(i2->first) != TGscat[i0->first][i1->first].end() &&
                     TGscat[i0->first][i1->first][i2->first].find(i4->first) != TGscat[i0->first][i1->first][i2->first].end() &&
                     TGscat[i0->first][i1->first][i2->first][i4->first].find(i5->first) != TGscat[i0->first][i1->first][i2->first][i4->first].end() &&
                     TGscat[i0->first][i1->first][i2->first][i4->first][i5->first].find(i6->first) != TGscat[i0->first][i1->first][i2->first][i4->first][i5->first].end() &&
                     TGscat[i0->first][i1->first][i2->first][i4->first][i5->first][i6->first].find(i7->first) != TGscat[i0->first][i1->first][i2->first][i4->first][i5->first][i6->first].end() &&
                     TGscat[i0->first][i1->first][i2->first][i4->first][i5->first][i6->first][i7->first].find(i8->first) != TGscat[i0->first][i1->first][i2->first][i4->first][i5->first][i6->first][i7->first].end()
                    )                 
                    {
//                     g = TGscat[i1->first][i2->first][i4->first];
                    }
                  else
                    {  
                     TGscat[i0->first][i1->first][i2->first][i4->first][i5->first][i6->first][i7->first][i8->first] = new TGraphErrors();
                    }  
                    g = TGscat[i0->first][i1->first][i2->first][i4->first][i5->first][i6->first][i7->first][i8->first]; 
                    g->SetPoint(g->GetN(), (i3->first)/10., (i8->second)/1. );
//                    g->Print(); 
                    g->SetLineWidth(3);
                   }//i8
                   }//i7
                  } //i6
//return 0;
                }  //i5
             }   //i4  
          } //i3
       } //i2
    } //i1
  } //i0

}

void book_graphics()
{
   leg = new TLegend(0.5,0.6,0.97,0.99); 
   leg->SetFillColor(0);

   lsmil = new TLegend(0.6,0.8,0.97,0.99); 
   leg->SetFillColor(0);


   iline02 = new TLine(2.4,1.,2.4,20.);
   iline02->SetLineColor(kRed);
   iline02->SetLineWidth(1);
   iline03 = new TLine(3.6,1.,3.6,20.);
   iline03->SetLineColor(kRed);
   iline03->SetLineWidth(1);
   iline05 = new TLine(6.0,1.,6.0,20.);
   iline05->SetLineColor(kRed);
   iline05->SetLineWidth(1);

   iline15 = new TLine(15.,1.,15.,20.);
   iline15->SetLineColor(kRed);
   iline15->SetLineWidth(1);

   iline120 = new TLine(120.,1.,120.,20.);
   iline120->SetLineColor(kRed);
   iline120->SetLineWidth(1);

   kline02 = new TLine(2.4,0.001,2.4,4.0);
   kline02->SetLineColor(kRed);
   kline02->SetLineWidth(1);
   kline03 = new TLine(3.6,0.001,3.6,4.0);
   kline03->SetLineColor(kRed);
   kline03->SetLineWidth(1);
   kline05 = new TLine(6.0,0.001,6.0,4.0);
   kline05->SetLineColor(kRed);
   kline05->SetLineWidth(1);

   kline15 = new TLine(15.,0.001,15.,4.0);
   kline15->SetLineColor(kRed);
   kline15->SetLineWidth(1);

   kline120 = new TLine(120.,0.001,120.,4.0);
   kline120->SetLineColor(kRed);
   kline120->SetLineWidth(1);


   gStyle->SetOptStat(0);
   gStyle->SetPadTopMargin(0.01);
   gStyle->SetPadRightMargin(0.03);
   gStyle->SetPadBottomMargin(0.1 );
   gStyle->SetPadLeftMargin(0.1);

   c1 = new TCanvas("c1","A Simple Graph with error bars",10,10,1000,1000);
   c1->SetGrid();

   h1=new TH1D("h1",";(beam momentum) p, GeV; #sigma_{pointing @ plane 3}, #mum",1,0.1,1000.);
   h1->SetMaximum( 50.);
   h1->SetMinimum(1.0);

   c2 = new TCanvas("c2","A Simple Graph with error bars",10,10,1000,1000);
   c2->SetGrid();

   h2 = new TH1D("h2",";(beam momentum) p, GeV; #sigma_{track slope @ plane 3}, #murad",1,0.1,1000.);
   h2->SetMaximum( 4.);
   h2->SetMinimum(0.001);

   c3 = new TCanvas("c3","A Simple Graph with error bars",10,10,1000,1000);
   c3->SetGrid();

   h3=new TH1D("h3",";z, mm; #sigma_{pointing @ plane 3}, #mum",1,-100.,1000.);
   h3->SetMaximum( 30.);
   h3->SetMinimum(0.0);

   c4 = new TCanvas("c4","A Simple Graph with error bars",10,10,1000,1000);
   c4->SetGrid();

   h4 = new TH1D("h4",";z, mm; #sigma_{pointing @ plane 3}, #mum",1,-20.,140.);
   h4->SetMaximum( 10.);
   h4->SetMinimum(0.00);

   c3->Divide(2,1);


   g = 0;
}

void plot_01(Int_t ithreshold, Int_t dsith, Int_t idkapton)

{
//   c1->GetFrame()->SetFillColor(21);
//   c1->GetFrame()->SetBorderSize(12);
//   c1->Divide(3,3);


   c1->cd();
   gPad->SetGrid();
   h1->DrawCopy();
   gPad->SetLogy(1);
   gPad->SetLogx(1);
   jc1=1;
   for(int i=20; i<=150;i+=130)
   {
    jc1--;
std::cout << " TGres[i][20][20][40][50][50][ipl][idkapton]; ipl="<<ipl<<" " <<  TGres[i][20][20][40][50][50][ipl][idkapton] << " " << std::endl;
    TGres[i][20][20][40][50][50][ipl][idkapton]->SetLineColor(100+jc1*10 );
    TGres[i][20][20][40][50][50][ipl][idkapton]->DrawClone("plsame");
   }

   jc1=1;
   for(int i=20; i<=150;i+=130)
   {
    jc1--;
    TGres[150][i][20][40][50][50][ipl][idkapton]->Print();
    TGres[150][i][20][40][50][50][ipl][idkapton]->SetLineColor( 80+jc1*10 );
    TGres[150][i][20][40][50][50][ipl][idkapton]->DrawClone("plsame");
   }

   jc1=1;
   for(int i=20; i<=150;i+=130)
   {
    jc1--;
    TGres[150][150][i][40][50][50][ipl][idkapton]->Print();
    TGres[150][150][i][40][50][50][ipl][idkapton]->SetLineColor( 60+jc1*10 );
    TGres[150][150][i][40][50][50][ipl][idkapton]->DrawClone("plsame");
   }
   iline02->Draw();iline03->Draw();iline05->Draw();
   iline15->Draw();
   iline120->Draw();
   leg->Clear();
   leg->AddEntry("","Seven telescope planes [0 to 6]:","");
   leg->AddEntry("","#Deltaz_{01} = #Deltaz_{56}","");
   leg->AddEntry("","#Deltaz_{12} = #Deltaz_{45}","");
   leg->AddEntry("","#Deltaz_{23} = #Deltaz_{34}","");
   leg->AddEntry("","configuration tags: #Deltaz_{01}#Deltaz_{12}#Deltaz_{23}","");
   leg->AddEntry("","with N [#Deltaz= 20 mm], W [#Deltaz= 150 mm]","");
   leg->AddEntry(TGres[20][20][20][40][50][50][ipl][idkapton],"NNN","l");
   leg->AddEntry(TGres[150][20][20][40][50][50][ipl][idkapton],"wNN","l");
   leg->AddEntry(TGres[150][150][20][40][50][50][ipl][idkapton],"wwN","l");
   leg->AddEntry(TGres[150][150][150][40][50][50][ipl][idkapton],"www","l");
   leg->Draw();
   c1->SaveAs("log_1_w_Material_presolution_01.pdf");
   c1->SaveAs("log_1_w_Material_presolution_01.png");
}

void plot_02(Int_t ithreshold, Int_t dsith, Int_t idkapton)
{
//   c1->cd(2);
   gPad->SetGrid();
//   h1->SetMaximum(10.);
   h1->DrawCopy();
   gPad->SetLogy(1);
   gPad->SetLogx(1);
   jc1=1;
   for(int i=20; i<=150;i+=130)
   {
    jc1--;
    TGres[20][i][20][40][50][50][ipl][idkapton]->SetLineColor(100 +jc1*10 );
    TGres[20][i][20][40][50][50][ipl][idkapton]->DrawClone("plsame");
   }

   jc1=1;
   for(int i=20; i<=150;i+=130)
   {
    jc1--;
    TGres[20][150][i][40][50][50][ipl][idkapton]->SetLineColor( 80 +jc1*10 );
    TGres[20][150][i][40][50][50][ipl][idkapton]->DrawClone("plsame");
   }
   iline02->Draw();iline03->Draw();iline05->Draw();
   iline15->Draw();
   iline120->Draw();
   leg->Clear();
   leg->AddEntry("","Seven telescope planes [0 to 6]:","");
   leg->AddEntry("","#Deltaz_{01} = #Deltaz_{56}","");
   leg->AddEntry("","#Deltaz_{12} = #Deltaz_{45}","");
   leg->AddEntry("","#Deltaz_{23} = #Deltaz_{34}","");
   leg->AddEntry("","configuration tags: #Deltaz_{01}#Deltaz_{12}#Deltaz_{23}","");
   leg->AddEntry("","with N [#Deltaz= 20 mm], W [#Deltaz= 150 mm]","");
   leg->AddEntry(TGres[20][20][20][40][50][50][ipl][idkapton],"NNN","l");
   leg->AddEntry(TGres[20][150][20][40][50][50][ipl][idkapton],"NwN","l");
   leg->AddEntry(TGres[20][150][150][40][50][50][ipl][idkapton],"Nww","l");
   leg->Draw();
   c1->SaveAs("log_1_w_Material_presolution_02.pdf");
   c1->SaveAs("log_1_w_Material_presolution_02.png");
}

void plot_03(Int_t ithreshold, Int_t dsith, Int_t idkapton)
{

//   c1->cd(3);
   gPad->SetGrid();
//   h1->SetMaximum(10.);
   h1->DrawCopy();
   gPad->SetLogy(1);
   gPad->SetLogx(1);
   jc1=1;
   for(int i=20; i<=150;i+=130)
   {
    jc1--;
    TGres[20][20][i][40][50][50][ipl][idkapton]->SetLineColor(100 +jc1*10  );
    TGres[20][20][i][40][50][50][ipl][idkapton]->DrawClone("plsame");
 
    TGres[20][i][150][40][50][50][ipl][idkapton]->SetLineColor(80 +jc1*10  );
    TGres[20][i][150][40][50][50][ipl][idkapton]->DrawClone("plsame");
 
    TGres[i][20][150][40][50][50][ipl][idkapton]->SetLineColor(60 +jc1*10  );
    TGres[i][20][150][40][50][50][ipl][idkapton]->DrawClone("plsame");
 
    TGres[i][150][150][40][50][50][ipl][idkapton]->SetLineColor(60 +jc1*10  );
    TGres[i][150][150][40][50][50][ipl][idkapton]->DrawClone("plsame");
  }
   iline02->Draw();iline03->Draw();iline05->Draw();
   iline15->Draw();
   iline120->Draw();
   leg->Clear();
   leg->AddEntry("","Seven telescope planes [0 to 6]:","");
   leg->AddEntry("","#Deltaz_{01} = #Deltaz_{56}","");
   leg->AddEntry("","#Deltaz_{12} = #Deltaz_{45}","");
   leg->AddEntry("","#Deltaz_{23} = #Deltaz_{34}","");
   leg->AddEntry("","configuration tags: #Deltaz_{01}#Deltaz_{12}#Deltaz_{23}","");
   leg->AddEntry("","with N [#Deltaz= 20 mm], W [#Deltaz= 150 mm]","");
   leg->AddEntry(TGres[20][20][20][40][50][50][ipl][idkapton],"NNN","l");
   leg->AddEntry(TGres[20][20][150][40][50][50][ipl][idkapton],"NNw","l");
   leg->AddEntry(TGres[20][150][150][40][50][50][ipl][idkapton],"Nww","l");
   leg->AddEntry(TGres[150][20][150][40][50][50][ipl][idkapton],"wNw","l");
   leg->AddEntry(TGres[150][150][150][40][50][50][ipl][idkapton],"www","l");
   leg->Draw();
   c1->SaveAs("log_1_w_Material_presolution_03.pdf");
   c1->SaveAs("log_1_w_Material_presolution_03.png");

}

void plot_04(Int_t ithreshold, Int_t dsith, Int_t idkapton)
{
  
//   c1->cd(4);
   gPad->SetGrid();
   h1->DrawCopy();
   gPad->SetLogy(1);
   gPad->SetLogx(1);
   jc1=1;
   for(int i=50; i<=50000;i=i*10)
   {  
    jc1--; 
    TGres[20][20][20][40][i][50][ipl][idkapton]->SetLineColor(kRed ); 
    TGres[20][20][20][40][i][50][ipl][idkapton]->SetLineStyle(1-jc1 ); 
    TGres[20][20][20][40][i][50][ipl][idkapton]->DrawClone("plsame");
    TGres[150][20][20][40][i][50][ipl][idkapton]->SetLineColor( kViolet );
    TGres[150][20][20][40][i][50][ipl][idkapton]->SetLineStyle(1-jc1 ); 
    TGres[150][20][20][40][i][50][ipl][idkapton]->DrawClone("plsame");
   }
   iline02->Draw();iline03->Draw();iline05->Draw();
   iline15->Draw();
   iline120->Draw();
   leg->Clear();
   leg->AddEntry("","Seven telescope planes [0 to 6]:","");
   leg->AddEntry("","#Deltaz_{01} = #Deltaz_{56}","");
   leg->AddEntry("","#Deltaz_{12} = #Deltaz_{45}","");
   leg->AddEntry("","#Deltaz_{23} = #Deltaz_{34}","");
   leg->AddEntry("","configuration tags: #Deltaz_{01}#Deltaz_{12}#Deltaz_{23}","");
   leg->AddEntry("","with N [#Deltaz= 20 mm], W [#Deltaz= 150 mm]","");
   leg->AddEntry(TGres[20][20][20][40][50][50][ipl][idkapton],"NNN DUT [Si] 50 #mum [thick]","l");
   leg->AddEntry(TGres[20][20][20][40][500][50][ipl][idkapton],"NNN DUT 500 #mum ","l");
   leg->AddEntry(TGres[20][20][20][40][5000][50][ipl][idkapton],"NNN DUT 5000 #mum","l");
   leg->AddEntry(TGres[20][20][20][40][50000][50][ipl][idkapton],"NNN DUT 50000 #mum","l");
   leg->AddEntry(TGres[150][20][20][40][50][50][ipl][idkapton],"wNN DUT 50 #mum","l");
   leg->AddEntry(TGres[150][20][20][40][500][50][ipl][idkapton],"wNN DUT 500 #mum","l");
   leg->AddEntry(TGres[150][20][20][40][5000][50][ipl][idkapton],"wNN DUT 5000 #mum","l");
   leg->AddEntry(TGres[150][20][20][40][50000][50][ipl][idkapton],"wNN DUT 50000 #mum","l");
   leg->Draw();
   c1->SaveAs("log_1_w_Material_presolution_04.pdf");
   c1->SaveAs("log_1_w_Material_presolution_04.png");

}

void plot_05(Int_t ithreshold, Int_t dsith, Int_t idkapton)
{

//   c1->cd(5);
   gPad->SetGrid();
   h1->DrawCopy();
   gPad->SetLogy(1);
   gPad->SetLogx(1);
   jc1=1;
   for(int i=50; i<=50000;i=i*10)
   {  
    jc1--; 
    TGres[20][20][20][40][i][50][ipl][idkapton]->SetLineColor(kRed ); 
    TGres[20][20][20][40][i][50][ipl][idkapton]->SetLineStyle(1-jc1 ); 
    TGres[20][20][20][40][i][50][ipl][idkapton]->DrawClone("plsame");
    TGres[20][150][20][40][i][50][ipl][idkapton]->SetLineColor( kGreen ) ; 
    TGres[20][150][20][40][i][50][ipl][idkapton]->SetLineStyle(1-jc1 ); 
    TGres[20][150][20][40][i][50][ipl][idkapton]->DrawClone("plsame");
   }
   iline02->Draw();iline03->Draw();iline05->Draw();
   iline15->Draw();
   iline120->Draw();
   leg->Clear();
   leg->AddEntry("","Seven telescope planes [0 to 6]:","");
   leg->AddEntry("","#Deltaz_{01} = #Deltaz_{56}","");
   leg->AddEntry("","#Deltaz_{12} = #Deltaz_{45}","");
   leg->AddEntry("","#Deltaz_{23} = #Deltaz_{34}","");
   leg->AddEntry("","configuration tags: #Deltaz_{01}#Deltaz_{12}#Deltaz_{23}","");
   leg->AddEntry("","with N [#Deltaz= 20 mm], W [#Deltaz= 150 mm]","");
   leg->AddEntry(TGres[20][20][20][40][50][50][ipl][idkapton],"NNN DUT [Si] 50 #mum [thick]","l");
   leg->AddEntry(TGres[20][20][20][40][500][50][ipl][idkapton],"NNN DUT 500 #mum ","l");
   leg->AddEntry(TGres[20][20][20][40][5000][50][ipl][idkapton],"NNN DUT 5000 #mum","l");
   leg->AddEntry(TGres[20][20][20][40][50000][50][ipl][idkapton],"NNN DUT 50000 #mum","l");
   leg->AddEntry(TGres[20][150][20][40][50][50][ipl][idkapton],"NwN DUT 50 #mum","l");
   leg->AddEntry(TGres[20][150][20][40][500][50][ipl][idkapton],"NwN DUT 500 #mum","l");
   leg->AddEntry(TGres[20][150][20][40][5000][50][ipl][idkapton],"NwN DUT 5000 #mum","l");
   leg->AddEntry(TGres[20][150][20][40][50000][50][ipl][idkapton],"NwN DUT 50000 #mum","l");
   leg->Draw();
   c1->SaveAs("log_1_w_Material_presolution_05.pdf");
   c1->SaveAs("log_1_w_Material_presolution_05.png");

}

void plot_06(Int_t ithreshold, Int_t dsith, Int_t idkapton)
{

//   c1->cd(6);
   gPad->SetGrid();
   h1->DrawCopy();
   gPad->SetLogy(1);
   gPad->SetLogx(1);
   jc1=1;
   for(int i=50; i<=50000;i=i*10)
   {  
    jc1--; 
    TGres[20][20][20][40][i][50][ipl][idkapton]->SetLineColor(kRed ); 
    TGres[20][20][20][40][i][50][ipl][idkapton]->SetLineStyle(1-jc1 ); 
    TGres[20][20][20][40][i][50][ipl][idkapton]->DrawClone("plsame");
    TGres[150][150][20][40][i][50][ipl][idkapton]->SetLineColor( kBlue ); 
    TGres[150][150][20][40][i][50][ipl][idkapton]->SetLineStyle(1-jc1 ); 
    TGres[150][150][20][40][i][50][ipl][idkapton]->DrawClone("plsame");
   }
   iline02->Draw();iline03->Draw();iline05->Draw();
   iline15->Draw();
   iline120->Draw();
   leg->Clear();
   leg->AddEntry("","Seven telescope planes [0 to 6]:","");
   leg->AddEntry("","#Deltaz_{01} = #Deltaz_{56}","");
   leg->AddEntry("","#Deltaz_{12} = #Deltaz_{45}","");
   leg->AddEntry("","#Deltaz_{23} = #Deltaz_{34}","");
   leg->AddEntry("","configuration tags: #Deltaz_{01}#Deltaz_{12}#Deltaz_{23}","");
   leg->AddEntry("","with N [#Deltaz= 20 mm], W [#Deltaz= 150 mm]","");
   leg->AddEntry(TGres[20][20][20][40][50][50][ipl][idkapton],"NNN DUT [Si] 50 #mum [thick]","l");
   leg->AddEntry(TGres[20][20][20][40][500][50][ipl][idkapton],"NNN DUT 500 #mum ","l");
   leg->AddEntry(TGres[20][20][20][40][5000][50][ipl][idkapton],"NNN DUT 5000 #mum","l");
   leg->AddEntry(TGres[20][20][20][40][50000][50][ipl][idkapton],"NNN DUT 50000 #mum","l");
   leg->AddEntry(TGres[150][150][20][40][50][50][ipl][idkapton],"wwN DUT 50 #mum","l");
   leg->AddEntry(TGres[150][150][20][40][500][50][ipl][idkapton],"wwN DUT 500 #mum","l");
   leg->AddEntry(TGres[150][150][20][40][5000][50][ipl][idkapton],"wwN DUT 5000 #mum","l");
   leg->AddEntry(TGres[150][150][20][40][50000][50][ipl][idkapton],"wwN DUT 50000 #mum","l");
   leg->Draw();
   c1->SaveAs("log_1_w_Material_presolution_06.pdf");
   c1->SaveAs("log_1_w_Material_presolution_06.png");
}

void plot_07(Int_t ithreshold, Int_t dsith, Int_t idkapton)
{
   gPad->SetGrid();
   h1->DrawCopy();
   gPad->SetLogy(1);
   gPad->SetLogx(1);
   jc1=1;
   for(int i=50; i<=50000;i=i*10)
   {  
    jc1--; 
    TGres[20][20][20][40][i][50][ipl][idkapton]->SetLineColor(kRed ); 
    TGres[20][20][20][40][i][50][ipl][idkapton]->SetLineStyle(1-jc1 ); 
    TGres[20][20][20][40][i][50][ipl][idkapton]->DrawClone("plsame");
    TGres[150][150][150][40][i][50][ipl][idkapton]->SetLineColor( kBlue ); 
    TGres[150][150][150][40][i][50][ipl][idkapton]->SetLineStyle(1-jc1 ); 
    TGres[150][150][150][40][i][50][ipl][idkapton]->DrawClone("plsame");
   }
   iline02->Draw();iline03->Draw();iline05->Draw();
   iline15->Draw();
   iline120->Draw();
   leg->Clear();
   leg->AddEntry("","Seven telescope planes [0 to 6]:","");
   leg->AddEntry("","#Deltaz_{01} = #Deltaz_{56}","");
   leg->AddEntry("","#Deltaz_{12} = #Deltaz_{45}","");
   leg->AddEntry("","#Deltaz_{23} = #Deltaz_{34}","");
   leg->AddEntry("","configuration tags: #Deltaz_{01}#Deltaz_{12}#Deltaz_{23}","");
   leg->AddEntry("","with N [#Deltaz= 20 mm], W [#Deltaz= 150 mm]","");
   leg->AddEntry(TGres[20][20][20][40][50][50][ipl][idkapton],"NNN DUT [Si] 50 #mum [thick]","l");
   leg->AddEntry(TGres[20][20][20][40][500][50][ipl][idkapton],"NNN DUT 500 #mum ","l");
   leg->AddEntry(TGres[20][20][20][40][5000][50][ipl][idkapton],"NNN DUT 5000 #mum","l");
   leg->AddEntry(TGres[20][20][20][40][50000][50][ipl][idkapton],"NNN DUT 50000 #mum","l");
   leg->AddEntry(TGres[150][150][150][40][50][50][ipl][idkapton],"www DUT 50 #mum","l");
   leg->AddEntry(TGres[150][150][150][40][500][50][ipl][idkapton],"www DUT 500 #mum","l");
   leg->AddEntry(TGres[150][150][150][40][5000][50][ipl][idkapton],"www DUT 5000 #mum","l");
   leg->AddEntry(TGres[150][150][150][40][50000][50][ipl][idkapton],"www DUT 50000 #mum","l");
   leg->Draw();
   c1->SaveAs("log_1_w_Material_presolution_07.pdf");
   c1->SaveAs("log_1_w_Material_presolution_07.png");

}


void plot_08(Int_t ithreshold, Int_t dsith, Int_t idkapton)
{

   gPad->SetGrid();
   h1->DrawCopy();
   gPad->SetLogy(1);
   gPad->SetLogx(1);
   jc1=2;
   for(int i=30; i<= 50;i+= 5)
   {
    jc1--;
    if( jc1 > 0 )
    {
      TGres[20][20][20][i][50][50][ipl][idkapton]->SetLineColor( 100  ); 
      TGres[150][150][150][i][50][50][ipl][idkapton]->SetLineColor( 50 ); 
      TGres[20][20][20][i][50][50][ipl][idkapton]->SetLineStyle( 2  ); 
      TGres[150][150][150][i][50][50][ipl][idkapton]->SetLineStyle( 2 ); 
    }
    else
    {
      TGres[20][20][20][i][50][50][ipl][idkapton]->SetLineColor( 100 +jc1*10   ); 
      TGres[150][150][150][i][50][50][ipl][idkapton]->SetLineColor( 50 +jc1*10 );  
      TGres[20][20][20][i][50][50][ipl][idkapton]->SetLineStyle( 1  ); 
      TGres[150][150][150][i][50][50][ipl][idkapton]->SetLineStyle( 1 ); 
    }
    TGres[20][20][20][i][50][50][ipl][idkapton]->DrawClone("plsame");
    TGres[150][150][150][i][50][50][ipl][idkapton]->DrawClone("plsame");
 
   }
   iline02->Draw();iline03->Draw();iline05->Draw();
   iline15->Draw();
   iline120->Draw();
   c1->SaveAs("log_1_w_Material_presolution_07.pdf");
   leg->Clear();
   leg->AddEntry("","Seven telescope planes [0 to 6]:","");
   leg->AddEntry("","#Deltaz_{01} = #Deltaz_{56}","");
   leg->AddEntry("","#Deltaz_{12} = #Deltaz_{45}","");
   leg->AddEntry("","#Deltaz_{23} = #Deltaz_{34}","");
   leg->AddEntry("","configuration tags: #Deltaz_{01}#Deltaz_{12}#Deltaz_{23}","");
   leg->AddEntry("","with N [#Deltaz= 20 mm], W [#Deltaz= 150 mm]","");
   leg->AddEntry(TGres[20][20][20][30][50][50][ipl][idkapton],"NNN, Mimosa26 #sigma_{intrinsic}=3.0 #mum ","l");
   leg->AddEntry(TGres[20][20][20][35][50][50][ipl][idkapton],"NNN, Mimosa26 #sigma_{intrinsic}=3.5 #mum ","l");
   leg->AddEntry(TGres[20][20][20][40][50][50][ipl][idkapton],"NNN, Mimosa26 #sigma_{intrinsic}=4.0 #mum ","l");
   leg->AddEntry(TGres[20][20][20][45][50][50][ipl][idkapton],"NNN, Mimosa26 #sigma_{intrinsic}=4.5 #mum ","l");
   leg->AddEntry(TGres[150][150][150][30][50][50][ipl][idkapton],"www, Mimosa26 #sigma_{intrinsic}=3.0 #mum","l");
   leg->AddEntry(TGres[150][150][150][35][50][50][ipl][idkapton],"www, Mimosa26 #sigma_{intrinsic}=3.5 #mum","l");
   leg->AddEntry(TGres[150][150][150][40][50][50][ipl][idkapton],"www, Mimosa26 #sigma_{intrinsic}=4.0 #mum","l");
   leg->AddEntry(TGres[150][150][150][45][50][50][ipl][idkapton],"www, Mimosa26 #sigma_{intrinsic}=4.5 #mum","l");
   leg->Draw();
   c1->SaveAs("log_1_w_Material_presolution_08.png");
   c1->SaveAs("log_1_w_Material_presolution_08.pdf");

}

void plot_09(Int_t ithreshold, Int_t dsith, Int_t idkapton)
{
//   c1->cd(9);
   gPad->SetGrid();
//   h1->SetMaximum(10.);
   h1->DrawCopy();
   gPad->SetLogy(1);
   gPad->SetLogx(1);
   jc1=1;
   for(int i=20; i<=150;i+=130)
   {
    for(int k=20; k<=150;k+=130)
    {
    jc1--;
    TGres[k][i][20][40][50][50][ipl][idkapton]->SetLineColor(100 +jc1*10  );
    TGres[k][i][20][40][50][50][ipl][idkapton]->DrawClone("plsame");

    TGres[k][i][150][40][50][50][ipl][idkapton]->SetLineColor(60 +jc1*10  );
    TGres[k][i][150][40][50][50][ipl][idkapton]->DrawClone("plsame");
    }
   }

   iline02->Draw();iline03->Draw();iline05->Draw();
   iline15->Draw();
   iline120->Draw();
   leg->Clear();
   leg->AddEntry("","Seven telescope planes [0 to 6]:","");
   leg->AddEntry("","#Deltaz_{01} = #Deltaz_{56}","");
   leg->AddEntry("","#Deltaz_{12} = #Deltaz_{45}","");
   leg->AddEntry("","#Deltaz_{23} = #Deltaz_{34}","");
   leg->AddEntry("","configuration tags: #Deltaz_{01}#Deltaz_{12}#Deltaz_{23}","");
   leg->AddEntry("","with N [#Deltaz= 20 mm], W [#Deltaz= 150 mm]","");
   leg->AddEntry(TGres[20][20][20][40][50][50][ipl][idkapton],"NNN ","l");
   leg->AddEntry(TGres[150][20][20][40][50][50][ipl][idkapton],"wNN ","l");
   leg->AddEntry(TGres[20][150][20][40][50][50][ipl][idkapton],"NwN ","l");
   leg->AddEntry(TGres[20][20][150][40][50][50][ipl][idkapton],"NNw ","l");
   leg->AddEntry(TGres[20][150][150][40][50][50][ipl][idkapton],"Nww ","l");
   leg->AddEntry(TGres[150][20][150][40][50][50][ipl][idkapton],"wNw ","l");
   leg->AddEntry(TGres[150][150][20][40][50][50][ipl][idkapton],"wwN ","l");
   leg->AddEntry(TGres[150][150][150][40][50][50][ipl][idkapton],"www ","l");
   leg->Draw();
   c1->SaveAs("log_1_w_Material_presolution_09.pdf");
   c1->SaveAs("log_1_w_Material_presolution_09.png");

}

void scat_01(Int_t ithreshold, Int_t dsith, Int_t idkapton)
{
   c2->cd();
   gPad->SetGrid();
   h2->DrawCopy();
   gPad->SetLogy(1);
   gPad->SetLogx(1);
   jc2=1;
   for(int i=20; i<=150;i+=130)
   {
    jc2--;
//       i0:20 i1 20 i2 20 i3 24 i4 50 i5 60 i6 60 i7 3 i8 2000 i8:sec: 0.118
std::cout << " " << i << " [20][20][40][50][50][ " << ipl << " ][ " << idkapton << std::endl;
    TGscat[i][20][20][40][50][50][ipl][idkapton]->Print();
    TGscat[i][20][20][40][50][50][ipl][idkapton]->SetLineColor(100+jc2*10 );
    TGscat[i][20][20][40][50][50][ipl][idkapton]->DrawClone("plsame");
   }

   jc2=1;
   for(int i=20; i<=150;i+=130)
   {
    jc2--;
    TGscat[150][i][20][40][50][50][ipl][idkapton]->SetLineColor( 80+jc2*10 );
    TGscat[150][i][20][40][50][50][ipl][idkapton]->DrawClone("plsame");
   }

   jc2=1;
   for(int i=20; i<=150;i+=130)
   {
    jc2--;
    TGscat[150][150][i][40][50][50][ipl][idkapton]->SetLineColor( 60+jc2*10 );
    TGscat[150][150][i][40][50][50][ipl][idkapton]->DrawClone("plsame");
   }
   kline02->Draw();kline03->Draw();kline05->Draw();
   kline15->Draw();
   kline120->Draw();
   leg->Clear();
   leg->AddEntry("","Seven telescope planes [0 to 6]:","");
   leg->AddEntry("","#Deltaz_{01} = #Deltaz_{56}","");
   leg->AddEntry("","#Deltaz_{12} = #Deltaz_{45}","");
   leg->AddEntry("","#Deltaz_{23} = #Deltaz_{34}","");
   leg->AddEntry("","configuration tags: #Deltaz_{01}#Deltaz_{12}#Deltaz_{23}","");
   leg->AddEntry("","with N [#Deltaz= 20 mm], W [#Deltaz= 150 mm]","");
   leg->AddEntry(TGscat[20][20][20][40][50][50][ipl][idkapton],"NNN","l");
   leg->AddEntry(TGscat[150][20][20][40][50][50][ipl][idkapton],"wNN","l");
   leg->AddEntry(TGscat[150][150][20][40][50][50][ipl][idkapton],"wwN","l");
   leg->AddEntry(TGscat[150][150][150][40][50][50][ipl][idkapton],"www","l");
   leg->Draw();
   c2->SaveAs("log_1_w_Material_scattering_01.pdf");
   c2->SaveAs("log_1_w_Material_scattering_01.png");

}

void scat_02(Int_t ithreshold, Int_t dsith, Int_t idkapton)
{
//   c2->cd(2);
   gPad->SetGrid();
//   h2->SetMaximum(10.);
   h2->DrawCopy();
   gPad->SetLogy(1);
   gPad->SetLogx(1);
   jc2=1;
   for(int i=20; i<=150;i+=130)
   {
    jc2--;
    TGscat[20][i][20][40][50][50][ipl][idkapton]->SetLineColor(100 +jc2*10 );
    TGscat[20][i][20][40][50][50][ipl][idkapton]->DrawClone("plsame");
   }

   jc2=1;
   for(int i=20; i<=150;i+=130)
   {
    jc2--;
    TGscat[20][150][i][40][50][50][ipl][idkapton]->SetLineColor( 80 +jc2*10 );
    TGscat[20][150][i][40][50][50][ipl][idkapton]->DrawClone("plsame");
   }
    kline02->Draw();kline03->Draw();kline05->Draw();
   kline15->Draw();
   kline120->Draw();
   leg->Clear();
   leg->AddEntry("","Seven telescope planes [0 to 6]:","");
   leg->AddEntry("","#Deltaz_{01} = #Deltaz_{56}","");
   leg->AddEntry("","#Deltaz_{12} = #Deltaz_{45}","");
   leg->AddEntry("","#Deltaz_{23} = #Deltaz_{34}","");
   leg->AddEntry("","configuration tags: #Deltaz_{01}#Deltaz_{12}#Deltaz_{23}","");
   leg->AddEntry("","with N [#Deltaz= 20 mm], W [#Deltaz= 150 mm]","");
   leg->AddEntry(TGscat[20][20][20][40][50][50][ipl][idkapton],"NNN","l");
   leg->AddEntry(TGscat[20][150][20][40][50][50][ipl][idkapton],"NwN","l");
   leg->AddEntry(TGscat[20][150][150][40][50][50][ipl][idkapton],"Nww","l");
   leg->Draw();
   c2->SaveAs("log_1_w_Material_scattering_02.pdf");
   c2->SaveAs("log_1_w_Material_scattering_02.png");

}

void scat_03(Int_t ithreshold, Int_t dsith, Int_t idkapton)
{

//   c2->cd(3);
   gPad->SetGrid();
//   h2->SetMaximum(10.);
   h2->DrawCopy();
   gPad->SetLogy(1);
   gPad->SetLogx(1);
   jc2=1;
   for(int i=20; i<=150;i+=130)
   {
    jc2--;
    TGscat[20][20][i][40][50][50][ipl][idkapton]->SetLineColor(100 +jc2*10  );
    TGscat[20][20][i][40][50][50][ipl][idkapton]->DrawClone("plsame");

    TGscat[20][i][150][40][50][50][ipl][idkapton]->SetLineColor( 80 +jc2*10  );
    TGscat[20][i][150][40][50][50][ipl][idkapton]->DrawClone("plsame");

    TGscat[i][20][150][40][50][50][ipl][idkapton]->SetLineColor( 60 +jc2*10  );
    TGscat[i][20][150][40][50][50][ipl][idkapton]->DrawClone("plsame");
   }
    kline02->Draw();kline03->Draw();kline05->Draw();
   kline15->Draw();
   kline120->Draw();
   leg->Clear();
   leg->AddEntry("","Seven telescope planes [0 to 6]:","");
   leg->AddEntry("","#Deltaz_{01} = #Deltaz_{56}","");
   leg->AddEntry("","#Deltaz_{12} = #Deltaz_{45}","");
   leg->AddEntry("","#Deltaz_{23} = #Deltaz_{34}","");
   leg->AddEntry("","configuration tags: #Deltaz_{01}#Deltaz_{12}#Deltaz_{23}","");
   leg->AddEntry("","with N [#Deltaz= 20 mm], W [#Deltaz= 150 mm]","");
   leg->AddEntry(TGscat[20][20][20][40][50][50][ipl][idkapton],"NNN","l");
   leg->AddEntry(TGscat[20][20][150][40][50][50][ipl][idkapton],"NNw","l");
   leg->AddEntry(TGscat[20][150][150][40][50][50][ipl][idkapton],"Nww","l");
   leg->AddEntry(TGscat[150][20][150][40][50][50][ipl][idkapton],"wNw","l");
   leg->Draw();
   c2->SaveAs("log_1_w_Material_scattering_03.pdf");
   c2->SaveAs("log_1_w_Material_scattering_03.png");


}

void scat_04(Int_t ithreshold, Int_t dsith, Int_t idkapton)
{
//   c2->cd(4);
   gPad->SetGrid();
   h2->DrawCopy();
   gPad->SetLogy(1);
   gPad->SetLogx(1);
   jc2=1;
   for(int i=50; i<=50000;i=i*10)
   {  
    jc2--; 
    TGscat[20][20][20][40][i][50][ipl][idkapton]->SetLineColor(kRed ); 
    TGscat[20][20][20][40][i][50][ipl][idkapton]->SetLineStyle(1-jc2 ); 
    TGscat[20][20][20][40][i][50][ipl][idkapton]->DrawClone("plsame");
    TGscat[150][20][20][40][i][50][ipl][idkapton]->SetLineColor( kViolet );
    TGscat[150][20][20][40][i][50][ipl][idkapton]->SetLineStyle(1-jc2 ); 
    TGscat[150][20][20][40][i][50][ipl][idkapton]->DrawClone("plsame");
   }
    kline02->Draw();kline03->Draw();kline05->Draw();
   kline15->Draw();
   kline120->Draw();
   leg->Clear();
   leg->AddEntry("","Seven telescope planes [0 to 6]:","");
   leg->AddEntry("","#Deltaz_{01} = #Deltaz_{56}","");
   leg->AddEntry("","#Deltaz_{12} = #Deltaz_{45}","");
   leg->AddEntry("","#Deltaz_{23} = #Deltaz_{34}","");
   leg->AddEntry("","configuration tags: #Deltaz_{01}#Deltaz_{12}#Deltaz_{23}","");
   leg->AddEntry("","with N [#Deltaz= 20 mm], W [#Deltaz= 150 mm]","");
   leg->AddEntry(TGscat[20][20][20][40][50][50][ipl][idkapton],"NNN DUT [Si] 50 #mum [thick]","l");
   leg->AddEntry(TGscat[20][20][20][40][500][50][ipl][idkapton],"NNN DUT 500 #mum ","l");
   leg->AddEntry(TGscat[20][20][20][40][5000][50][ipl][idkapton],"NNN DUT 5000 #mum","l");
   leg->AddEntry(TGscat[20][20][20][40][50000][50][ipl][idkapton],"NNN DUT 50000 #mum","l");
   leg->AddEntry(TGscat[150][20][20][40][50][50][ipl][idkapton],"wNN DUT 50 #mum","l");
   leg->AddEntry(TGscat[150][20][20][40][500][50][ipl][idkapton],"wNN DUT 500 #mum","l");
   leg->AddEntry(TGscat[150][20][20][40][5000][50][ipl][idkapton],"wNN DUT 5000 #mum","l");
   leg->AddEntry(TGscat[150][20][20][40][50000][50][ipl][idkapton],"wNN DUT 50000 #mum","l");
   leg->Draw();
   c2->SaveAs("log_1_w_Material_scattering_04.pdf");
   c2->SaveAs("log_1_w_Material_scattering_04.png");

}

void scat_05(Int_t ithreshold, Int_t dsith, Int_t idkapton)
{
//   c2->cd(5);
   gPad->SetGrid();
   h2->DrawCopy();
   gPad->SetLogy(1);
   gPad->SetLogx(1);
   jc2=1;
   for(int i=50; i<=50000;i=i*10)
   {  
    jc2--; 
    TGscat[20][20][20][40][i][50][ipl][idkapton]->SetLineColor(kRed ); 
    TGscat[20][20][20][40][i][50][ipl][idkapton]->SetLineStyle(1-jc2 ); 
    TGscat[20][20][20][40][i][50][ipl][idkapton]->DrawClone("plsame");
    TGscat[20][150][20][40][i][50][ipl][idkapton]->SetLineColor( kGreen ); 
    TGscat[20][150][20][40][i][50][ipl][idkapton]->SetLineStyle(1-jc2 ); 
    TGscat[20][150][20][40][i][50][ipl][idkapton]->DrawClone("plsame"); 
 
   }
   kline02->Draw();kline03->Draw();kline05->Draw();
   kline15->Draw();
   kline120->Draw();
   leg->Clear();
   leg->AddEntry("","Seven telescope planes [0 to 6]:","");
   leg->AddEntry("","#Deltaz_{01} = #Deltaz_{56}","");
   leg->AddEntry("","#Deltaz_{12} = #Deltaz_{45}","");
   leg->AddEntry("","#Deltaz_{23} = #Deltaz_{34}","");
   leg->AddEntry("","configuration tags: #Deltaz_{01}#Deltaz_{12}#Deltaz_{23}","");
   leg->AddEntry("","with N [#Deltaz= 20 mm], W [#Deltaz= 150 mm]","");
   leg->AddEntry(TGscat[20][20][20][40][50][50][ipl][idkapton],"NNN DUT [Si] 50 #mum [thick]","l");
   leg->AddEntry(TGscat[20][20][20][40][500][50][ipl][idkapton],"NNN DUT 500 #mum ","l");
   leg->AddEntry(TGscat[20][20][20][40][5000][50][ipl][idkapton],"NNN DUT 5000 #mum","l");
   leg->AddEntry(TGscat[20][20][20][40][50000][50][ipl][idkapton],"NNN DUT 50000 #mum","l");
   leg->AddEntry(TGscat[20][150][20][40][50][50][ipl][idkapton],"NwN DUT 50 #mum","l");
   leg->AddEntry(TGscat[20][150][20][40][500][50][ipl][idkapton],"NwN DUT 500 #mum","l");
   leg->AddEntry(TGscat[20][150][20][40][5000][50][ipl][idkapton],"NwN DUT 5000 #mum","l");
   leg->AddEntry(TGscat[20][150][20][40][50000][50][ipl][idkapton],"NwN DUT 50000 #mum","l");
   leg->Draw();
   c2->SaveAs("log_1_w_Material_scattering_05.pdf");
   c2->SaveAs("log_1_w_Material_scattering_05.png");

}

void scat_06(Int_t ithreshold, Int_t dsith, Int_t idkapton)
{
//   c2->cd(6);
   gPad->SetGrid();
   h2->DrawCopy();
   gPad->SetLogy(1);
   gPad->SetLogx(1);
   jc2=1;
   for(int i=50; i<=50000;i=i*10)
   {  
    jc2--; 
    TGscat[20][20][20][40][i][50][ipl][idkapton]->SetLineColor(kRed ); 
    TGscat[20][20][20][40][i][50][ipl][idkapton]->SetLineStyle(1 -jc2 ); 
    TGscat[20][20][20][40][i][50][ipl][idkapton]->DrawClone("plsame");
    TGscat[150][150][20][40][i][50][ipl][idkapton]->SetLineColor( kBlue ); 
    TGscat[150][150][20][40][i][50][ipl][idkapton]->SetLineStyle(1 -jc2 ); 
    TGscat[150][150][20][40][i][50][ipl][idkapton]->DrawClone("plsame");
   }
   kline02->Draw();kline03->Draw();kline05->Draw();
   kline15->Draw();
   kline120->Draw();
   leg->Clear();
   leg->AddEntry("","Seven telescope planes [0 to 6]:","");
   leg->AddEntry("","#Deltaz_{01} = #Deltaz_{56}","");
   leg->AddEntry("","#Deltaz_{12} = #Deltaz_{45}","");
   leg->AddEntry("","#Deltaz_{23} = #Deltaz_{34}","");
   leg->AddEntry("","configuration tags: #Deltaz_{01}#Deltaz_{12}#Deltaz_{23}","");
   leg->AddEntry("","with N [#Deltaz= 20 mm], W [#Deltaz= 150 mm]","");
   leg->AddEntry(TGscat[20][20][20][40][50][50][ipl][idkapton],"NNN DUT [Si] 50 #mum [thick]","l");
   leg->AddEntry(TGscat[20][20][20][40][500][50][ipl][idkapton],"NNN DUT 500 #mum ","l");
   leg->AddEntry(TGscat[20][20][20][40][5000][50][ipl][idkapton],"NNN DUT 5000 #mum","l");
   leg->AddEntry(TGscat[20][20][20][40][50000][50][ipl][idkapton],"NNN DUT 50000 #mum","l");
   leg->AddEntry(TGscat[150][150][20][40][50][50][ipl][idkapton],"wwN DUT 50 #mum","l");
   leg->AddEntry(TGscat[150][150][20][40][500][50][ipl][idkapton],"wwN DUT 500 #mum","l");
   leg->AddEntry(TGscat[150][150][20][40][5000][50][ipl][idkapton],"wwN DUT 5000 #mum","l");
   leg->AddEntry(TGscat[150][150][20][40][50000][50][ipl][idkapton],"wwN DUT 50000 #mum","l");
   leg->Draw();
   c2->SaveAs("log_1_w_Material_scattering_06.pdf");
   c2->SaveAs("log_1_w_Material_scattering_06.png");


}

void scat_07(Int_t ithreshold, Int_t dsith, Int_t idkapton)
{
// c2->cd(7);
   gPad->SetGrid();
   h2->DrawCopy();
   gPad->SetLogy(1);
   gPad->SetLogx(1);
   jc2=1;
   for(int i=50; i<=50000;i=i*10)
   {  
    jc2--; 
    TGscat[20][20][20][40][i][50][ipl][idkapton]->SetLineColor(kRed ); 
    TGscat[20][20][20][40][i][50][ipl][idkapton]->SetLineStyle(1 -jc2 ); 
    TGscat[20][20][20][40][i][50][ipl][idkapton]->DrawClone("plsame");
    TGscat[150][150][150][40][i][50][ipl][idkapton]->SetLineColor( kBlue ); 
    TGscat[150][150][150][40][i][50][ipl][idkapton]->SetLineStyle(1 -jc2 ); 
    TGscat[150][150][150][40][i][50][ipl][idkapton]->DrawClone("plsame");
   }
//   gPad->Print("log_1_w_Material_scattering_9.png");

   kline02->Draw();kline03->Draw();kline05->Draw();
   kline15->Draw();
   kline120->Draw();
   leg->Clear();
   leg->AddEntry("","Seven telescope planes [0 to 6]:","");
   leg->AddEntry("","#Deltaz_{01} = #Deltaz_{56}","");
   leg->AddEntry("","#Deltaz_{12} = #Deltaz_{45}","");
   leg->AddEntry("","#Deltaz_{23} = #Deltaz_{34}","");
   leg->AddEntry("","configuration tags: #Deltaz_{01}#Deltaz_{12}#Deltaz_{23}","");
   leg->AddEntry("","with N [#Deltaz= 20 mm], W [#Deltaz= 150 mm]","");
   leg->AddEntry(TGscat[20][20][20][40][50][50][ipl][idkapton],"NNN DUT [Si] 50 #mum [thick]","l");
   leg->AddEntry(TGscat[20][20][20][40][500][50][ipl][idkapton],"NNN DUT 500 #mum ","l");
   leg->AddEntry(TGscat[20][20][20][40][5000][50][ipl][idkapton],"NNN DUT 5000 #mum","l");
   leg->AddEntry(TGscat[20][20][20][40][50000][50][ipl][idkapton],"NNN DUT 50000 #mum","l");
   leg->AddEntry(TGscat[150][150][150][40][50][50][ipl][idkapton],"www DUT 50 #mum","l");
   leg->AddEntry(TGscat[150][150][150][40][500][50][ipl][idkapton],"www DUT 500 #mum","l");
   leg->AddEntry(TGscat[150][150][150][40][5000][50][ipl][idkapton],"www DUT 5000 #mum","l");
   leg->AddEntry(TGscat[150][150][150][40][50000][50][ipl][idkapton],"www DUT 50000 #mum","l");
   leg->Draw();
   c2->SaveAs("log_1_w_Material_scattering_07.pdf");
   c2->SaveAs("log_1_w_Material_scattering_07.png");

}

void scat_08(Int_t ithreshold, Int_t dsith, Int_t idkapton)
{
//   c2->cd(8);
   gPad->SetGrid();
   h2->DrawCopy();
   gPad->SetLogy(1);
   gPad->SetLogx(1);
   jc2=2;
   for(int i=30; i<= 50;i+= 5)
   {
    jc2--;
    if( jc2 > 0 )
    {
      TGscat[20][20][20][i][50][50][ipl][idkapton]->SetLineColor( 100  ); 
      TGscat[150][150][150][i][50][50][ipl][idkapton]->SetLineColor( 50 ); 
      TGscat[20][20][20][i][50][50][ipl][idkapton]->SetLineStyle( 2  ); 
      TGscat[150][150][150][i][50][50][ipl][idkapton]->SetLineStyle( 2 ); 
    }
    else
    {
      TGscat[20][20][20][i][50][50][ipl][idkapton]->SetLineColor( 100 +jc2*10   ); 
      TGscat[150][150][150][i][50][50][ipl][idkapton]->SetLineColor( 50 +jc2*10 );  
      TGscat[20][20][20][i][50][50][ipl][idkapton]->SetLineStyle( 1  ); 
      TGscat[150][150][150][i][50][50][ipl][idkapton]->SetLineStyle( 1 ); 
    }
    TGscat[20][20][20][i][50][50][ipl][idkapton]->DrawClone("plsame");
    TGscat[150][150][150][i][50][50][ipl][idkapton]->DrawClone("plsame");
 
   }
   kline02->Draw();kline03->Draw();kline05->Draw();
   kline15->Draw();
   kline120->Draw();
   leg->Clear();
   leg->AddEntry("","Seven telescope planes [0 to 6]:","");
   leg->AddEntry("","#Deltaz_{01} = #Deltaz_{56}","");
   leg->AddEntry("","#Deltaz_{12} = #Deltaz_{45}","");
   leg->AddEntry("","#Deltaz_{23} = #Deltaz_{34}","");
   leg->AddEntry("","configuration tags: #Deltaz_{01}#Deltaz_{12}#Deltaz_{23}","");
   leg->AddEntry("","with N [#Deltaz= 20 mm], W [#Deltaz= 150 mm]","");
   leg->AddEntry(TGscat[20][20][20][30][50][50][ipl][idkapton],"NNN, Mimosa26 #sigma_{intrinsic}=3.0 #mum ","l");
   leg->AddEntry(TGscat[20][20][20][35][50][50][ipl][idkapton],"NNN, Mimosa26 #sigma_{intrinsic}=3.5 #mum ","l");
   leg->AddEntry(TGscat[20][20][20][40][50][50][ipl][idkapton],"NNN, Mimosa26 #sigma_{intrinsic}=4.0 #mum ","l");
   leg->AddEntry(TGscat[20][20][20][45][50][50][ipl][idkapton],"NNN, Mimosa26 #sigma_{intrinsic}=4.5 #mum ","l");
   leg->AddEntry(TGscat[150][150][150][30][50][50][ipl][idkapton],"www, Mimosa26 #sigma_{intrinsic}=3.0 #mum","l");
   leg->AddEntry(TGscat[150][150][150][35][50][50][ipl][idkapton],"www, Mimosa26 #sigma_{intrinsic}=3.5 #mum","l");
   leg->AddEntry(TGscat[150][150][150][40][50][50][ipl][idkapton],"www, Mimosa26 #sigma_{intrinsic}=4.0 #mum","l");
   leg->AddEntry(TGscat[150][150][150][45][50][50][ipl][idkapton],"www, Mimosa26 #sigma_{intrinsic}=4.5 #mum","l");
   leg->Draw();
   c2->SaveAs("log_1_w_Material_scattering_08.pdf");
   c2->SaveAs("log_1_w_Material_scattering_08.png");

}

void scat_09(Int_t ithreshold, Int_t dsith, Int_t idkapton)
{
//   c2->cd(9);
   gPad->SetGrid();
//   h2->SetMaximum(10.);
   h2->DrawCopy();
   gPad->SetLogy(1);
   gPad->SetLogx(1);
   jc2=1;
   for(int i=20; i<=150;i+=130)
   {
    for(int k=20; k<=150;k+=130)
    {
    jc2--;
      TGscat[k][i][20][40][50][50][ipl][idkapton]->SetLineColor(100 +jc2*10  );
      TGscat[k][i][20][40][50][50][ipl][idkapton]->DrawClone("plsame"); 
 
      TGscat[k][i][150][40][50][50][ipl][idkapton]->SetLineColor(60 +jc2*10  );
      TGscat[k][i][150][40][50][50][ipl][idkapton]->DrawClone("plsame"); 
    }
   }
   kline02->Draw();kline03->Draw();kline05->Draw();
   kline15->Draw();
   kline120->Draw();
   leg->Clear();
   leg->AddEntry("","Seven telescope planes [0 to 6]:","");
   leg->AddEntry("","#Deltaz_{01} = #Deltaz_{56}","");
   leg->AddEntry("","#Deltaz_{12} = #Deltaz_{45}","");
   leg->AddEntry("","#Deltaz_{23} = #Deltaz_{34}","");
   leg->AddEntry("","configuration tags: #Deltaz_{01}#Deltaz_{12}#Deltaz_{23}","");
   leg->AddEntry("","with N [#Deltaz= 20 mm], W [#Deltaz= 150 mm]","");
   leg->AddEntry(TGscat[20][20][20][40][50][50][ipl][idkapton],"NNN","l");
   leg->AddEntry(TGscat[150][20][20][40][50][50][ipl][idkapton],"wNN","l");
   leg->AddEntry(TGscat[20][150][20][40][50][50][ipl][idkapton],"NwN","l");
   leg->AddEntry(TGscat[20][20][150][40][50][50][ipl][idkapton],"NNw","l");
   leg->AddEntry(TGscat[20][150][150][40][50][50][ipl][idkapton],"Nww","l");
   leg->AddEntry(TGscat[150][150][20][40][50][50][ipl][idkapton],"wwN","l");
   leg->AddEntry(TGscat[150][20][150][40][50][50][ipl][idkapton],"wNw","l");
   leg->AddEntry(TGscat[150][150][150][40][50][50][ipl][idkapton],"www","l");
   leg->Draw();
   c2->SaveAs("log_1_w_Material_scattering_09.pdf");
   c2->SaveAs("log_1_w_Material_scattering_09.png");

}

void reset_histograms()
{
 for(int isensor=0; isensor<6;isensor+=1)
 for(int ithickness=50; ithickness<=50000;ithickness+=10)
 for(int isith=10; isith<=200;ithickness+=10)
 for(int i1=20; i1<=150;i1+=130)
 for(int i2=20; i2<=150;i2+=130)
 for(int i3=20; i3<=150;i3+=130)
 for(int i4=30; i4<= 50;i4+= 1)
 for(int idkapton=250; idkapton<= 500;idkapton *= 2)
   {
    TGres[i1][i2][i3][i4][ithickness][isith][isensor][idkapton]->SetLineStyle(1); 
   }  

}

void plot_smilie_pbeam()
{

// global loop over DATA thresholds 3 to 12

for(int ithr=5;ithr<13;ithr++)
{
std::cout << "ithr" << ithr << std::endl;

c1->cd();


for(int resthr=35;resthr<45;resthr=resthr+1)
{
int dsith=50;
int idkapton=250;

//   c1->cd(10);
   gPad->SetGrid();
//   h1->SetMaximum(10.);
   h1->DrawCopy();
   gPad->SetLogy(1);
//   gPad->SetLogy(0);
   gPad->SetLogx(1);

int iline = 0;
int iplane=3;

   std::cout << "iplane" << iplane << std::endl;

   iline02->Draw();iline03->Draw();iline05->Draw();
   iline15->Draw();
   iline120->Draw();
   leg->Clear();
   leg->AddEntry("","Seven telescope planes [0 to 6]:","");
   leg->AddEntry("","#Deltaz_{01} = #Deltaz_{56}","");
   leg->AddEntry("","#Deltaz_{12} = #Deltaz_{45}","");
   leg->AddEntry("","#Deltaz_{23} = #Deltaz_{34}","");
   leg->AddEntry("","configuration tags: #Deltaz_{01}#Deltaz_{12}#Deltaz_{23}","");
   leg->AddEntry("","with N [#Deltaz= 20 mm], W [#Deltaz= 150 mm]","");

   iline++;


   jc2=0;
 
   int  iloop[] = {10,20,40,80,150};
   std::vector<int> vloop (iloop, iloop + sizeof(iloop) / sizeof(iloop[0]) );   

for(int idist=10;idist<160;idist=idist+10)
{
   bool skip= true;
   if( std::find( vloop.begin(), vloop.end(), idist ) != vloop.end() ) skip=false;
   if(skip)continue;

   std::cout << " " << idist  << " " <<  jc2 << std::endl;
   jc2++;


   TGres[idist][idist][idist][resthr][dsith][dsith][iplane][idkapton]->SetLineColor(jc2 );
   TGres[idist][idist][idist][resthr][dsith][dsith][iplane][idkapton]->SetLineStyle( 1 );
   TGres[idist][idist][idist][resthr][dsith][dsith][iplane][idkapton]->DrawClone("plsame");

   leg->AddEntry(TGres[idist][idist][idist] [resthr][dsith][dsith][iplane][idkapton],"--- ","l");
}


   jc1=1;
   for( t_dz_beam_iter i0  = measurement.begin();
                      i0 != measurement.end();
                      i0++)
   { 
    jc1--; 
    for( t_beam_thr_iter i1  = i0->second.begin();
                    i1 != i0->second.end();
                    i1++)
     {
       jc2=0; 
       for( t_thr_plane_iter i2  = i1->second.begin();
                      i2 != i1->second.end();
                      i2++)
       {
        if( i2->first != ithr  ) continue; 
        for( t_plane_meas_iter i3  = i2->second.begin();
                      i3 != i2->second.end();
                      i3++)
          {
            if( i3->first != iplane  ) continue; 
//            std::cout  << " " << i0->first << " " << i2->first<< " " << i3->first << std::endl;
            TGmeas[i0->first][i2->first][i3->first]->SetMarkerStyle(20+jc2++);        
            TGmeas[i0->first][i2->first][i3->first]->SetMarkerColor(kRed+jc1);        
            TGmeas[i0->first][i2->first][i3->first]->SetMarkerSize(1);        
            TGmeas[i0->first][i2->first][i3->first]->Draw("P same"); 
          }
       }
     }
   }



   leg->Draw();
   TString s1="log_1_w_Material_presolution_thr_"+SSTR(ithr)+"_"+SSTR(resthr);//+"_"+SSTR(iplane);
   c1->SaveAs(s1+".pdf");
   c1->SaveAs(s1+".png");
}
}

}


void plot_smilie_150(Int_t ithreshold, Int_t dsith, Int_t idkapton)
{

for(int resthr=30;resthr<=50;resthr=resthr+1)
{

c3->cd(1);
h3->DrawCopy();
c3->cd(2);
h4->DrawCopy();


for(int ithr=ithreshold;ithr<=ithreshold ;ithr++)
{
   std::cout << "ithr" << ithr << std::endl;
//   c1->cd(10);
   gPad->SetGrid();
//   h1->SetMaximum(10.);
   gPad->SetLogy(0);
   gPad->SetLogx(0);

   int iline = 0;

   jc2=0;
for(int i=0;i<rbeam.size();i++)
{
   jc2++;
   lsmil->Clear();

   iline++;
  
   if( !(rbeam[i]==24 || rbeam[i]==36 || rbeam[i] ==60) )continue;
   c3->cd(1);
//
   std::cout << " TSres [150][150][150]["<<resthr<<"]["<<dsith<<"]["<<dsith<<"]["<<rbeam[i]<<"]["<<idkapton<<"] " 
   <<   TSres[150][150][150][resthr][dsith][dsith][rbeam[i]][idkapton] << std::endl;
   TSres[150][150][150][resthr][dsith][dsith][rbeam[i]][idkapton]->SetLineColor(100  );
   TSres[150][150][150][resthr][dsith][dsith][rbeam[i]][idkapton]->SetLineStyle(iline   );
   TSres[150][150][150][resthr][dsith][dsith][rbeam[i]][idkapton]->DrawClone("plsame");

   for(int jc3=0;jc3<vbeam.size();jc3++)
   {
    if(vbeam[jc3] != rbeam[i])continue;
    h= TSmeas[150][ithr][vbeam[jc3]];
    h->SetMarkerStyle(20+jc2++);        
    h->SetMarkerColor(kRed+jc1);        
    h->SetMarkerSize(1);        
    h->Draw("P same"); 
    lsmil->AddEntry(h,"NNN ","l");
   }
//
   c3->cd(2);
//
   TSres[ 20][ 20][ 20][resthr][dsith][dsith][rbeam[i]][idkapton]->SetLineColor(100  );
   TSres[ 20][ 20][ 20][resthr][dsith][dsith][rbeam[i]][idkapton]->SetLineStyle(iline   );
   TSres[ 20][ 20][ 20][resthr][dsith][dsith][rbeam[i]][idkapton]->DrawClone("plsame");

   for(int jc3=0;jc3<vbeam.size();jc3++)
   {
    if(vbeam[jc3] != rbeam[i])continue;
    h= TSmeas[ 20][ithr][vbeam[jc3]];
    h->SetMarkerStyle(20+jc2++);        
    h->SetMarkerColor(kRed+jc1);        
    h->SetMarkerSize(1);        
    h->Draw("P same"); 
    lsmil->AddEntry(h,"NNN ","l");
   }


}
   lsmil->Draw();
   TString s3="smilie_xxx_dkapton_"+SSTR(idkapton)+"_sith_"+SSTR(dsith)+"_res_"+SSTR(resthr)+"_thr_"+SSTR(ithr);//+"_"+SSTR(iplane);
   c3->SaveAs(s3+".pdf");
   c3->SaveAs(s3+".png");
}
}

}
/*
void plot_smilie_20(Int_t ithreshold, Int_t dsith)
{
c4->cd();

for(int resthr=30;resthr<=49;resthr=resthr+1)
{
h4->DrawCopy();

for(int ithr=ithreshold;ithr<=ithreshold ;ithr++)
{
std::cout << "ithr" << ithr << std::endl;
//   c1->cd(10);
   gPad->SetGrid();
//   h1->SetMaximum(10.);
   gPad->SetLogy(0);
   gPad->SetLogx(0);

int iline = 0;

   jc2=0;
for(int i=0;i<rbeam.size();i++)
{
//std::cout << "i: " << rbeam[i] << std::endl;
   lsmil->Clear();
   if( !(rbeam[i]==24 || rbeam[i]==36 || rbeam[i] ==60) )continue;
 
   iline++;

   TSres[20][20][20][resthr][dsith][dsith][rbeam[i]][idkapton]->SetLineColor(100  );
   TSres[20][20][20][resthr][dsith][dsith][rbeam[i]][idkapton]->SetLineStyle(iline   );
   TSres[20][20][20][resthr][dsith][dsith][rbeam[i]][idkapton]->DrawClone("plsame");

   for(int jc3=0;jc3<vbeam.size();jc3++)
   {
    if(vbeam[jc3] != rbeam[i])continue;
    h= TSmeas[20][ithr][vbeam[jc3]];
//              h->Print(); 
//     std::cout << "20 i: " << rbeam[i] << " .vs. " << vbeam[jc3] << std::endl;
             
    h->SetMarkerStyle(20+jc2++);        
    h->SetMarkerColor(kRed+jc1);        
    h->SetMarkerSize(1);        
    h->Draw("P same"); 
    lsmil->AddEntry(h,"NNN ","l");
   }

}
   lsmil->Draw();
   TString s1="smilie_20_presolution_"+SSTR(dsith)+"_res_"+SSTR(resthr)+"_thr_"+SSTR(ithr);//+"_"+SSTR(iplane);
   c4->SaveAs(s1+".pdf");
   c4->SaveAs(s1+".png");


}
}

}

*/
int main(int argc, const char* argv)
{

read_DATURA();
read_gbl();
book_graphics();

// for ipl

// ipl=3; //inner plane GLOBAL

plot_smilie_pbeam();

for(int i0= 4;i0<=12;i0=i0+1)
{
for(int i1=50;i1<= 50;i1=i1+10)
{
for(int i2=250;i2<= 250;i2*=2 )
{
// plot_smilie_150(i0,i1,i2);
// plot_smilie_20(i0,i1);

 ipl= 3;
/*
 plot_01(i0,i1,i2);
 plot_02(i0,i1,i2);
 plot_03(i0,i1,i2);
 plot_04(i0,i1,i2);
 plot_05(i0,i1,i2);
 plot_06(i0,i1,i2);
 plot_07(i0,i1,i2); 
 plot_08(i0,i1,i2);
 plot_09(i0,i1,i2);
*/
// scattering:
// scat_01(i0,i1,i2);
// scat_02(i0,i1,i2);
// scat_03(i0,i1,i2);
/* scat_04(i0,i1,i2);
 scat_05(i0,i1,i2);
 scat_06(i0,i1,i2);
 scat_07(i0,i1,i2);*/
// scat_08(i0,i1,i2);
// scat_09(i0,i1,i2);

}
}
}




return 0;
}



