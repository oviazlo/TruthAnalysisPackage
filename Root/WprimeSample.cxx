// WprimeSample.cxx

#include <TruthAnalysis/WprimeSample.h>

#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

WprimeSample::WprimeSample(Name salias)
: m_alias(salias), m_name(salias),
  m_genxspb(0), m_calxspb(0), m_run(0) {
  // If needed, conver alias to a differnt name here.
  // Fill remaining properties based on assigned name.
  // W BG's.
  if ( m_name == "walpgnp0" ){
    m_title = "W #rightarrow e#nu + 0 jets";
    m_genxspb = 6932;
    m_calxspb = 8284;
    m_latex = "\\wenu";
    m_run = 107680;
  } else if ( m_name == "walpgnp1" ){
    m_title = "W #rightarrow e#nu + 1 jets";
    m_genxspb = 1306;
    m_calxspb = 1562;
    m_latex = "\\wenu";
    m_run = 107681;
  } else if ( m_name == "walpgnp2" ){
    m_title = "W #rightarrow e#nu + 2 jets";
    m_genxspb = 378;
    m_calxspb = 453;
    m_latex = "\\wenu";
    m_run = 107682;
  } else if ( m_name == "walpgnp3" ){
    m_title = "W #rightarrow e#nu + 3 jets";
    m_genxspb = 102;
    m_calxspb = 122;
    m_latex = "\\wenu";
    m_run = 107683;
  } else if ( m_name == "walpgnp4" ){
    m_title = "W #rightarrow e#nu + 4 jets";
    m_genxspb = 25.7;
    m_calxspb = 31.0;
    m_latex = "\\wenu";
    m_run = 107684;
  } else if ( m_name == "walpgnp5" ){
    m_title = "W #rightarrow e#nu + 5 jets";
    m_genxspb = 7.0;
    m_calxspb = 8.3;
    m_latex = "\\wenu";
    m_run = 107685;
  } else if ( m_name == "wenu" ) {
     m_title = "W #rightarrow e #nu";
     m_genxspb =  8938;
     m_calxspb = 10460;
     m_latex = "\\wenu";
     m_run = 106043;
  } else if ( m_name == "wplusenu" ) {
    m_title = "W^{+} #rightarrow e #nu";
    m_genxspb = 6891.0;
    m_calxspb = 7124.0;
    m_latex = "\\wplusenu";
    m_run = 147800;
 } else if ( m_name == "wminenu" ) {
    m_title = "W^{-} #rightarrow e #nu";
    m_genxspb = 4790.2;
    m_calxspb = 5048.0;
    m_latex = "\\wminenu";
    m_run = 147803;
  } else if ( m_name == "wmunu" ) {
     m_title = "W #rightarrow #mu #nu";
     m_genxspb =  8938;
     m_calxspb = 10460;
     m_latex = "\\wmunu";
     m_run = 106044;
  } else if (m_name == "wplusmunu"){
    m_title = "W^{+} #rightarrow #mu #nu";
    m_genxspb = 6891.0;
    m_calxspb = 7124.0;
    m_latex = "\\wplusmunu";
    m_run = 147801;
  } else if (m_name == "wminmunu"){
    m_title = "W^{-} #rightarrow #mu #nu";
    m_genxspb = 4790.2;
    m_calxspb = 5048.0;
    m_latex = "\\wminmunu";
    m_run = 147804;
  } else if ( m_name == "wtaunu" ) {
    m_title = "W #rightarrow #tau #nu";
    m_genxspb =  8699.0;
    m_calxspb =  10180.0;
    m_latex = "\\wtaunu";
    m_run = 107054;
  } else if (m_name == "wplustaunu"){
    m_title = "W^{+} #rightarrow #tau #nu";
    m_genxspb = 6889.8;
    m_calxspb = 6889.8;
    m_latex = "\\wplustaunu";
    m_run = 147802;
  } else if (m_name == "wmintaunu"){
    m_title = "W^{-} #rightarrow #tau #nu";
    m_genxspb = 4790.9;
    m_calxspb = 4790.9;
    m_latex = "\\wmintaunu";
    m_run = 147805;
  } else if ( m_name == "wlnu_200_500" ) {
     m_title = "W #rightarrow l #nu (200 < m < 500 GeV)";
     m_genxspb =  8.90;
     m_calxspb = 10.42;
     m_latex = "\\wrlnu[200,500]";
     m_run = 106604;
  } else if ( m_name == "wlnu_500_1500" ) {
     m_title = "W #rightarrow l #nu (500 < m < 1500 GeV)";
     m_genxspb = 0.242;
     m_calxspb = 0.283;
     m_latex = "\\wrlnu[500,1500]";
     m_run = 105274;
  } else if ( m_name == "wlnu_1500_2500" ) {
     m_title = "W #rightarrow l #nu (1500 < m < 2500 GeV)";
     m_genxspb = 0.00079;
     m_calxspb = 0.00092;
     m_latex = "\\wrlnu[1500,2500]";
     m_run = 105275;
  } else if ( m_name == "wlnu_2500" ) {
     m_title = "W #rightarrow l #nu (m > 2500 GeV)";
     m_genxspb = 0.00000904;
     m_calxspb = 0.0000072;  // Rough estimate: K=0.8 from fig 99, 2010 backup
     m_latex = "\\wrlnu[2500+]";
     m_run = 105276;
  } else if ( m_name == "wplusenu_200_500" ) {
    m_title = "W^{+} #rightarrow e #nu (200 < m < 500 GeV)";
    m_genxspb = 2.507;
    m_calxspb = 2.507;
    m_latex = "\\wplusrenu[200,500]";
    m_run = 129561;
  } else if ( m_name == "wminenu_200_500" ) {
    m_title = "W^{-} #rightarrow e #nu (200 < m < 500 GeV)";
    m_genxspb = 1.448;
    m_calxspb = 1.448;
    m_latex = "\\wminrenu[200,500]";
    m_run = 129566;
  } else if ( m_name == "wplusmunu_200_500" ) {
    m_title = "W^{+} #rightarrow #mu #nu (200 < m < 500 GeV)";
    m_genxspb = 2.507;
    m_calxspb = 2.507;
    m_latex = "\\wplusrmunu[200,500]";
    m_run = 129571;
  } else if ( m_name == "wminmunu_200_500" ) {
    m_title = "W^{-} #rightarrow #mu #nu (200 < m < 500 GeV)";
    m_genxspb = 1.448;
    m_calxspb = 1.448;
    m_latex = "\\wminrmunu[200,500]";
    m_run = 129576;
  } else if ( m_name == "wplustaunu_200_500" ) {
    m_title = "W^{+} #rightarrow #tau #nu (200 < m < 500 GeV)";
    m_genxspb = 2.507;
    m_calxspb = 2.507;
    m_latex = "\\wplusrtaunu[200,500]";
    m_run = 129581;
  } else if ( m_name == "wmintaunu_200_500" ) {
    m_title = "W^{-} #rightarrow #tau #nu (200 < m < 500 GeV)";
    m_genxspb = 1.448;
    m_calxspb = 1.448;
    m_latex = "\\wminrtaunu[200,500]";
    m_run = 129586;
  } else if ( m_name == "wplusenu_500_1500" ) {
    m_title = "W^{+} #rightarrow e #nu (500 < m < 1500 GeV)";
    m_genxspb = 0.07672;
    m_calxspb = 0.07672;
    m_latex = "\\wplusrenu[500,1500]";
    m_run = 129562;
  } else if ( m_name == "wminenu_500_1500" ) {
    m_title = "W^{-} #rightarrow e #nu (500 < m < 1500 GeV)";
    m_genxspb = 0.03423;
    m_calxspb = 0.03423;
    m_latex = "\\wminrenu[500,1500]";
    m_run = 129567;
  } else if ( m_name == "wplusmunu_500_1500" ) {
    m_title = "W^{+} #rightarrow #mu #nu (500 < m < 1500 GeV)";
    m_genxspb = 0.07672;
    m_calxspb = 0.07672;
    m_latex = "\\wplusrmunu[500,1500]";
    m_run = 129572;
  } else if ( m_name == "wminmunu_500_1500" ) {
    m_title = "W^{-} #rightarrow #mu #nu (500 < m < 1500 GeV)";
    m_genxspb = 0.03423;
    m_calxspb = 0.03423;
    m_latex = "\\wminrmunu[500,1500]";
    m_run = 129577;
  } else if ( m_name == "wplustaunu_500_1500" ) {
    m_title = "W^{+} #rightarrow #tau #nu (500 < m < 1500 GeV)";
    m_genxspb = 0.07672;
    m_calxspb = 0.07672;
    m_latex = "\\wplusrtaunu[500,1500]";
    m_run = 129582;
  } else if ( m_name == "wmintaunu_500_1500" ) {
    m_title = "W^{-} #rightarrow #tau #nu (500 < m < 1500 GeV)";
    m_genxspb = 0.03423;
    m_calxspb = 0.03423;
    m_latex = "\\wminrtaunu[500,1500]";
    m_run = 129587;
  } else if ( m_name == "wplusenu_1500_2500" ) {
    m_title = "W^{+} #rightarrow e #nu (1500 < m < 2500 GeV)";
    m_genxspb = 0.0003196;
    m_calxspb = 0.0003196;
    m_latex = "\\wplusrenu[1500,2500]";
    m_run = 129563;
  } else if ( m_name == "wminenu_1500_2500" ) {
    m_title = "W^{-} #rightarrow e #nu (1500 < m < 2500 GeV)";
    m_genxspb = 0.00009978;
    m_calxspb = 0.00009978;
    m_latex = "\\wminrenu[1500,2500]";
    m_run = 129568;
  } else if ( m_name == "wplusmunu_1500_2500" ) {
    m_title = "W^{+} #rightarrow #mu #nu (1500 < m < 2500 GeV)";
    m_genxspb = 0.0003196;
    m_calxspb = 0.0003196;
    m_latex = "\\wplusrmunu[1500,2500]";
    m_run = 129573;
  } else if ( m_name == "wminmunu_1500_2500" ) {
    m_title = "W^{-} #rightarrow #mu #nu (1500 < m < 2500 GeV)";
    m_genxspb = 0.00009978;
    m_calxspb = 0.00009978;
    m_latex = "\\wminrmunu[1500,2500]";
    m_run = 129578;
  } else if ( m_name == "wplustaunu_1500_2500" ) {
    m_title = "W^{+} #rightarrow #tau #nu (1500 < m < 2500 GeV)";
    m_genxspb = 0.0003196;
    m_calxspb = 0.0003196;
    m_latex = "\\wplusrtaunu[1500,2500]";
    m_run = 129583;
  } else if ( m_name == "wmintaunu_1500_2500" ) {
    m_title = "W^{-} #rightarrow #tau #nu (1500 < m < 2500 GeV)";
    m_genxspb = 0.00009978;
    m_calxspb = 0.00009978;
    m_latex = "\\wminrtaunu[1500,2500]";
    m_run = 129588;
  } else if ( m_name == "wplusenu_2500_3000" ) {
    m_title = "W^{+} #rightarrow e #nu (2500 < m < 3000 GeV)";
    m_genxspb = 0.000004369;
    m_calxspb = 0.000004369;
    m_latex = "\\wplusrenu[2500,3000]";
    m_run = 129564;
  } else if ( m_name == "wminenu_2500_3000" ) {
    m_title = "W^{-} #rightarrow e #nu (2500 < m < 3000 GeV)";
    m_genxspb = 0.000001465;
    m_calxspb = 0.000001465;
    m_latex = "\\wminrenu[2500,3000]";
    m_run = 129569;
  } else if ( m_name == "wplusmunu_2500_3000" ) {
    m_title = "W^{+} #rightarrow #mu #nu (2500 < m < 3000 GeV)";
    m_genxspb = 0.000004369;
    m_calxspb = 0.000004369;
    m_latex = "\\wplusrmunu[2500,3000]";
    m_run = 129574;
  } else if ( m_name == "wminmunu_2500_3000" ) {
    m_title = "W^{-} #rightarrow #mu #nu (2500 < m < 3000 GeV)";
    m_genxspb = 0.000001465;
    m_calxspb = 0.000001465;
    m_latex = "\\wminrmunu[2500,3000]";
    m_run = 129579;
  } else if ( m_name == "wplustaunu_2500_3000" ) {
    m_title = "W^{+} #rightarrow #tau #nu (2500 < m < 3000 GeV)";
    m_genxspb = 0.000004369;
    m_calxspb = 0.000004369;
    m_latex = "\\wplusrtaunu[2500,3000]";
    m_run = 129584;
  } else if ( m_name == "wmintaunu_2500_3000" ) {
    m_title = "W^{-} #rightarrow #tau #nu (2500 < m < 3000 GeV)";
    m_genxspb = 0.000001465;
    m_calxspb = 0.000001465;
    m_latex = "\\wminrtaunu[2500,3000]";
    m_run = 129589;
  } else if ( m_name == "wplusenu_3000_3500" ) {
    m_title = "W^{+} #rightarrow e #nu (3000 < m < 3500 GeV)";
    m_genxspb = 0.0000005541;
    m_calxspb = 0.0000005541;
    m_latex = "\\wplusrenu[3000,3500]";
    m_run = 129600;
  } else if ( m_name == "wminenu_3000_3500" ) {
    m_title = "W^{-} #rightarrow e #nu (3000 < m < 3500 GeV)";
    m_genxspb = 0.0000002132;
    m_calxspb = 0.0000002132;
    m_latex = "\\wminrenu[3000,3500]";
    m_run = 129602;
  } else if ( m_name == "wplusmunu_3000_3500" ) {
    m_title = "W^{+} #rightarrow #mu #nu (3000 < m < 3500 GeV)";
    m_genxspb = 0.0000005541;
    m_calxspb = 0.0000005541;
    m_latex = "\\wplusrmunu[3000,3500]";
    m_run = 129604;
  } else if ( m_name == "wminmunu_3000_3500" ) {
    m_title = "W^{-} #rightarrow #mu #nu (3000 < m < 3500 GeV)";
    m_genxspb = 0.0000002132;
    m_calxspb = 0.0000002132;
    m_latex = "\\wminrmunu[3000,3500]";
    m_run = 129606;
  } else if ( m_name == "wplustaunu_3000_3500" ) {
    m_title = "W^{+} #rightarrow #tau #nu (3000 < m < 3500 GeV)";
    m_genxspb = 0.0000005541;
    m_calxspb = 0.0000005541;
    m_latex = "\\wplusrtaunu[3000,3500]";
    m_run = 129608;
  } else if ( m_name == "wmintaunu_3000_3500" ) {
    m_title = "W^{-} #rightarrow #tau #nu (3000 < m < 3500 GeV)";
    m_genxspb = 0.0000002132;
    m_calxspb = 0.0000002132;
    m_latex = "\\wminrtaunu[3000,3500]";
    m_run = 129610;
  } else if ( m_name == "wplusenu_3500" ) {
    m_title = "W^{+} #rightarrow e #nu (m > 3500 GeV)";
    m_genxspb = 0.0000000774;
    m_calxspb = 0.0000000774;
    m_latex = "\\wplusrenu[3500]";
    m_run = 129601;
  } else if ( m_name == "wminenu_3500" ) {
    m_title = "W^{-} #rightarrow e #nu (m > 3500 GeV)";
    m_genxspb = 0.00000003461;
    m_calxspb = 0.00000003461;
    m_latex = "\\wminrenu[3500]";
    m_run = 129603;
  } else if ( m_name == "wplusmunu_3500" ) {
    m_title = "W^{+} #rightarrow #mu #nu (m > 3500 GeV)";
    m_genxspb = 0.0000000774;
    m_calxspb = 0.0000000774;
    m_latex = "\\wplusrmunu[3500]";
    m_run = 129605;
  } else if ( m_name == "wminmunu_3500" ) {
    m_title = "W^{-} #rightarrow #mu #nu (m > 3500 GeV)";
    m_genxspb = 0.00000003461;
    m_calxspb = 0.00000003461;
    m_latex = "\\wminrmunu[3500]";
    m_run = 129607;
  } else if ( m_name == "wplustaunu_3500" ) {
    m_title = "W^{+} #rightarrow #tau #nu (m > 3500 GeV)";
    m_genxspb = 0.0000000774;
    m_calxspb = 0.0000000774;
    m_latex = "\\wplusrtaunu[3500]";
    m_run = 129609;
  } else if ( m_name == "wmintaunu_3500" ) {
    m_title = "W^{-} #rightarrow #tau #nu (m > 3500 GeV)";
    m_genxspb = 0.00000003461;
    m_calxspb = 0.00000003461;
    m_latex = "\\wminrtaunu[3500]";
    m_run = 129611;
  // W Lepton Filtered ========
  } else if ( m_name == "wminenu_pt100" ) {
    m_title = "W^{-} #rightarrow e #nu (p_{T}^{e} > 100 GeV)";
    m_latex = "\\wminenu\\ ($\\pte > 100 \\gev$)";
    m_genxspb = 13.243; //64.130 * 0.20651
    m_calxspb = 13.243;
    m_run = 129593;
  } else if ( m_name == "wplusenu_pt100" ) {
    m_title = "W^{+} #rightarrow e #nu (p_{T}^{e} > 100 GeV)";
    m_latex = "\\wplusenu\\ ($\\pte > 100 \\gev$)";
    m_genxspb = 13.842; //69.270 * 0.19982
    m_calxspb = 13.842;
    m_run = 129590;
  } else if ( m_name == "wminenu_pt150" ) {
    m_title = "W^{-} #rightarrow e #nu (p_{T}^{e} > 150 GeV)";
    m_latex = "\\wminenu\\ ($\\pte > 150 \\gev$)";
    m_genxspb = 2.075; //64.130 * 0.032360
    m_calxspb = 2.075;
    m_run = 129623;
  } else if ( m_name == "wplusenu_pt150" ) {
    m_title = "W^{+} #rightarrow e #nu (p_{T}^{e} > 150 GeV)";
    m_latex = "\\wplusenu\\ ($\\pte > 150 \\gev$)";
    m_genxspb = 2.164; //69.270 * 0.031246
    m_calxspb = 2.164;
    m_run = 129620;
  } else if ( m_name == "wminenu_pt200" ) {
    m_title = "W^{-} #rightarrow e #nu (p_{T}^{e} > 200 GeV)";
    m_latex = "\\wminenu\\ ($\\pte > 200 \\gev$)";
    m_genxspb = 0.6390; // 2.5820 * 0.24749
    m_calxspb = 0.6390;
    m_run = 129633;
  } else if ( m_name == "wplusenu_pt200" ) {
    m_title = "W^{+} #rightarrow e #nu (p_{T}^{e} > 200 GeV)";
    m_latex = "\\wplusenu\\ ($\\pte > 200 \\gev$)";
    m_genxspb = 0.6674; // 2.6620 * 0.25073
    m_calxspb = 0.6674;
    m_run = 129630;
  } else if ( m_name == "wminenu_pt300" ) {
    m_title = "W^{-} #rightarrow e #nu (p_{T}^{e} > 300 GeV)";
    m_latex = "\\wminenu\\ ($\\pte > 300 \\gev$)";
    m_genxspb = 0.09954; // 2.5820 * 0.038551
    m_calxspb = 0.09954;
    m_run = 129643;
  } else if ( m_name == "wplusenu_pt300" ) {
    m_title = "W^{+} #rightarrow e #nu (p_{T}^{e} > 300 GeV)";
    m_latex = "\\wplusenu\\ ($\\pte > 300 \\gev$)";
    m_genxspb = 0.09954; // 2.5820 * 0.038551
    m_calxspb = 0.09954;
    m_run = 129640;
  } else if ( m_name == "wminpt100_mu" ) {
    m_title = "W^{-} #rightarrow #mu #nu (p_{T}^{#mu} > 100 GeV)";
    m_latex = "\\wminmunu\\ ($\\ptmu > 100 \\gev$)";
    m_genxspb = 13.719; //64.130 * 0.21393
    m_calxspb = 13.719;
    m_run = 129594;
  } else if ( m_name == "wpluspt100_mu" ) {
    m_title = "W^{+} #rightarrow #mu #nu (p_{T}^{#mu} > 100 GeV)";
    m_latex = "\\wplusmunu\\ ($\\ptmu > 100 \\gev$)";
    m_genxspb = 14.394; //69.270 * 0.20779
    m_calxspb = 14.394;
    m_run = 129591;
  } else if ( m_name == "wminpt150_mu" ) {
    m_title = "W^{-} #rightarrow #mu #nu (p_{T}^{#mu} > 150 GeV)";
    m_latex = "\\wminmunu\\ ($\\ptmu > 150 \\gev$)";
    m_genxspb = 2.1663; //64.130 * 0.033780
    m_calxspb = 2.1663;
    m_run = 129624;
  } else if ( m_name == "wpluspt150_mu" ) {
    m_title = "W^{+} #rightarrow #mu #nu (p_{T}^{#mu} > 150 GeV)";
    m_latex = "\\wplusmunu\\ ($\\pte > 150 \\gev$)";
    m_genxspb = 2.2579; //69.270 * 0.032596
    m_calxspb = 2.2579;
    m_run = 129621;
  } else if ( m_name == "wminpt200_mu" ) {
    m_title = "W^{-} #rightarrow #mu #nu (p_{T}^{#mu} > 200 GeV)";
    m_latex = "\\wminmunu\\ ($\\ptmu > 200 \\gev$)";
    m_genxspb = 0.66817; // 2.5820 * 0.25878
    m_calxspb = 0.66817;
    m_run = 129634;
  } else if ( m_name == "wpluspt200_mu" ) {
    m_title = "W^{+} #rightarrow #mu #nu (p_{T}^{#mu} > 200 GeV)";
    m_latex = "\\wplusmunu\\ ($\\ptmu > 200 \\gev$)";
    m_genxspb = 0.69499; // 2.6620 * 0.26108
    m_calxspb = 0.69499;
    m_run = 129631;
  } else if ( m_name == "wminpt300_mu" ) {
    m_title = "W^{-} #rightarrow #mu #nu (p_{T}^{#mu} > 300 GeV)";
    m_latex = "\\wminmunu\\ ($\\ptmu > 300 \\gev$)";
    m_genxspb = 0.10502; // 2.5820 * 0.040673
    m_calxspb = 0.10502;
    m_run = 129644;
  } else if ( m_name == "wpluspt300_mu" ) {
    m_title = "W^{+} #rightarrow #mu #nu (p_{T}^{#mu} > 300 GeV)";
    m_latex = "\\wplusenu\\ ($\\pte > 300 \\gev$)";
    m_genxspb = 0.11681; // 2.6620 * 0.043879
    m_calxspb = 0.11681;
    m_run = 129641;


  // Z BG's.===========================
  } else if ( m_name == "ztautau" ) {
    m_title = "Z #rightarrow #tau #tau";    
    m_genxspb = 1109.9;
    m_calxspb = 1109.9;
    m_latex = "\\ztautau";
    m_run = 147808;
  } else if ( m_name == "zee" ) {
    m_title = "Z #rightarrow e e";
    m_genxspb = 1109.9;
    m_calxspb = 1109.9;
    m_latex = "\\zee";
    m_run = 147806;
  } else if ( m_name == "zee_250_400" ) {
     m_title = "Z #rightarrow e e (250 < m < 400 GeV)";
     m_genxspb = 0.5492;
     m_calxspb = 0.5492;
     m_latex = "\\zree[250,400]";
     m_run = 129506;
  } else if ( m_name == "zee_400_600" ) {
     m_title = "Z #rightarrow e e (400 < m < 600 GeV)";
     m_latex = "\\zree[400,600]";
     m_genxspb = 0.08966;
     m_calxspb = 0.08966;
     m_run = 129507;
  } else if ( m_name == "zee_600_800" ) {
     m_title = "Z #rightarrow e e (600 < m < 800 GeV)";
     m_latex = "\\zree[600,800]";
     m_genxspb = 0.0151;
     m_calxspb = 0.0151;
     m_run = 129508;
  } else if ( m_name == "zee_800_1000" ) {
     m_title = "Z #rightarrow e e (800 < m < 1000 GeV)";
     m_latex = "\\zree[800,1000]";
     m_genxspb = 0.00375;
     m_calxspb = 0.00375;
     m_run = 129509;
  } else if ( m_name == "zee_1000_1250" ) {
     m_title = "Z #rightarrow e e (1000 < m < 1250 GeV)";
     m_latex = "\\zree[1000,1250]";
     m_genxspb = 0.001293;
     m_calxspb = 0.001293;
     m_run = 129510;
  } else if ( m_name == "zee_1250_1500" ) {
     m_title = "Z #rightarrow e e (1250 < m < 1500 GeV)";
     m_latex = "\\zree[1250,1500]";
     m_genxspb = 0.0003577;
     m_calxspb = 0.0003577;
     m_run = 129511;
  } else if ( m_name == "zee_1500_1750" ) {
     m_title = "Z #rightarrow e e (1500 < m < 1750 GeV)";
     m_latex = "\\zree[1500,1750]";
     m_genxspb = 0.0001123;
     m_calxspb = 0.0001123;
     m_run = 129512;
  } else if ( m_name == "zee_1750_2000" ) {
     m_title = "Z #rightarrow e e (1750 < m < 2000 GeV)";
     m_latex = "\\zree[1750,2000]";
     m_genxspb = 0.000038380;
     m_calxspb = 0.000038380;
     m_run = 129513;
  } else if ( m_name == "zee_2000_2250" ) {
     m_title = "Z #rightarrow e e (2000 < m < 2250 GeV)";
     m_latex = "\\zree[2000,2250]";
     m_genxspb = 0.000013890;
     m_calxspb = 0.000013890;
     m_run = 129514;
  } else if ( m_name == "zee_2250_2500" ) {
     m_title = "Z #rightarrow e e (2250 < m < 2500 GeV)";
     m_latex = "\\zree[2250,2500]";
     m_genxspb = 0.0000052260;
     m_calxspb = 0.0000052260;
     m_run = 129515;
  } else if ( m_name == "zee_2500_2750" ) {
     m_title = "Z #rightarrow e e (2500 < m < 2750 GeV)";
     m_latex = "\\zree[2500,2750]";
     m_genxspb = 0.0000020170;
     m_calxspb = 0.0000020170;
     m_run = 129516;
  } else if ( m_name == "zee_2750_3000" ) {
     m_title = "Z #rightarrow e e (2750 < m < 3000 GeV)";
     m_latex = "\\zree[2750,3000]";
     m_genxspb = 0.00000078910;
     m_calxspb = 0.00000078910;
     m_run = 129517;
  } else if ( m_name == "zee_3000" ) {
     m_title = "Z #rightarrow e e (m > 3000 GeV)";
     m_latex = "\\zree[3000+]";
     m_genxspb = 0.00000050390;
     m_calxspb = 0.00000050390;
     m_run = 129518;
  } else if ( m_name == "zmumu" ) {
     m_title = "Z #rightarrow #mu #mu";
     m_latex = "\\zmumu";
     m_genxspb = 1109.8;
     m_calxspb = 1109.8;
     m_run = 147807;
  } else if ( m_name == "zmumu_250_400" ) {
     m_title = "Z #rightarrow #mu #mu (250 < m < 400 GeV)";
     m_latex = "\\zrmumu[250,400]";
     m_genxspb = 0.54920;
     m_calxspb = 0.54920;
     m_run = 129526;
  } else if ( m_name == "zmumu_400_600" ) {
     m_title = "Z #rightarrow #mu #mu (400 < m < 600 GeV)";
     m_latex = "\\zrmumu[400,600]";
     m_genxspb = 0.08966;
     m_calxspb = 0.08966;
     m_run = 129527;
  } else if ( m_name == "zmumu_600_800" ) {
     m_title = "Z #rightarrow #mu #mu (600 < m < 800 GeV)";
     m_latex = "\\zrmumu[600,800]";
     m_genxspb = 0.015100;
     m_calxspb = 0.015100;
     m_run = 129528;
  } else if ( m_name == "zmumu_800_1000" ) {
     m_title = "Z #rightarrow #mu #mu (800 < m < 1000 GeV)";
     m_latex = "\\zrmumu[800,1000]";
     m_genxspb = 0.0037500;
     m_calxspb = 0.0037500;
     m_run = 129529;
  } else if ( m_name == "zmumu_1000_1250" ) {
     m_title = "Z #rightarrow #mu #mu (1000 < m < 1250 GeV)";
     m_latex = "\\zrmumu[1000,1250]";
     m_genxspb = 0.0012930;
     m_calxspb = 0.0012930;
     m_run = 129530;
  } else if ( m_name == "zmumu_1250_1500" ) {
     m_title = "Z #rightarrow #mu #mu (1250 < m < 1500 GeV)";
     m_latex = "\\zrmumu[1250,1500]";
     m_genxspb = 0.0003577;
     m_calxspb = 0.0003577;
     m_run = 129531;
  } else if ( m_name == "zmumu_1500_1750" ) {
     m_title = "Z #rightarrow #mu #mu (1500 < m < 1750 GeV)";
     m_latex = "\\zrmumu[1500,1750]";
     m_genxspb = 0.00011230;
     m_calxspb = 0.00011230;
     m_run = 129532;
  } else if ( m_name == "zmumu_1750_2000" ) {
     m_title = "Z #rightarrow #mu #mu (1750 < m < 2000 GeV)";
     m_latex = "\\zrmumu[1750,2000]";
     m_genxspb = 0.000038380;
     m_calxspb = 0.000038380;
     m_run = 129533;
  } else if ( m_name == "zmumu_2000_2250" ) {
     m_title = "Z #rightarrow #mu #mu (2000 < m < 2250 GeV)";
     m_latex = "\\zrmumu[2000,2250]";
     m_genxspb = 0.00001389;
     m_calxspb = 0.00001389;
     m_run = 129534;
  } else if ( m_name == "zmumu_2250_2500" ) {
    m_title = "Z #rightarrow #mu #mu (2250 < m < 2500 GeV)";
    m_latex = "\\zrmumu[2250,2500]";
    m_genxspb = 0.0000052260; 
    m_calxspb = 0.0000052260; 
    m_run = 129535;
  } else if ( m_name == "zmumu_2500_2750" ) {
    m_title = "Z #rightarrow #mu #mu (2500 < m < 2750 GeV)";
    m_latex = "\\zrmumu[2500,2750]";
    m_genxspb = 0.0000020170;
    m_calxspb = 0.0000020170; 
    m_run = 129536;
  } else if ( m_name == "zmumu_2750_3000" ) {
    m_title = "Z #rightarrow #mu #mu (2750 < m < 3000 GeV)";
    m_latex = "\\zrmumu[2750,3000]";
    m_genxspb = 0.00000078910; 
    m_calxspb = 0.00000078910; 
    m_run = 129537;
  } else if ( m_name == "zmumu_3000" ) {
    m_title = "Z #rightarrow #mu #mu (m > 3000 GeV)";
    m_latex = "\\zrmumu[3000+]";
    m_genxspb = 0.0000005039; // 
    m_calxspb = 0.0000005039; // 
    m_run = 129538;

  } else if ( m_name == "ztautau_250_400" ) {
    m_title = "Z #rightarrow #tau #tau (250 < m < 400 GeV)";
    m_genxspb = 0.5492;
    m_calxspb = 0.5492;
    m_latex = "\\zrtautau[250,400]";
    m_run = 129546;
  } else if ( m_name == "ztautau_400_600" ) {
    m_title = "Z #rightarrow #tau #tau (400 < m < 600 GeV)";
    m_genxspb = 0.08966;
    m_calxspb = 0.08966;
    m_latex = "\\zrtautau[400,600]";
    m_run = 129547;
  } else if ( m_name == "ztautau_600_800" ) {
    m_title = "Z #rightarrow #tau #tau (600 < m <800 GeV)";
    m_genxspb = 0.0151;
    m_calxspb = 0.0151;
    m_latex = "\\zrtautau[600,800]";
    m_run = 129548;
  } else if ( m_name == "ztautau_800_1000" ) {
    m_title = "Z #rightarrow #tau #tau (800 < m < 1000 GeV)";
    m_genxspb = 0.00375;
    m_calxspb = 0.00375;
    m_latex = "\\zrtautau[800,1000]";
    m_run = 129549;
  } else if ( m_name == "ztautau_1000_1250" ) {
    m_title = "Z #rightarrow #tau #tau (1000 < m < 1250 GeV)";
    m_genxspb = 0.001293;
    m_calxspb = 0.001293;
    m_latex = "\\zrtautau[1000,1250]";
    m_run = 129550;
  } else if ( m_name == "ztautau_1250_1500" ) {
    m_title = "Z #rightarrow #tau #tau (1250 < m < 1500 GeV)";
    m_genxspb = 0.0003577;
    m_calxspb = 0.0003577;
    m_latex = "\\zrtautau[1250,1500]";
    m_run = 129551;
  } else if (m_name == "ztautau_1500_1750" ) {
    m_title = "Z #rightarrow #tau #tau (1500 < m < 1750 GeV)";
    m_genxspb = 0.0001123;
    m_calxspb = 0.0001123;
    m_latex = "\\zrtautau[1500,1750]";
    m_run = 129552;
  } else if ( m_name == "ztautau_1750_2000" ) {
    m_title = "Z #rightarrow #tau #tau (1750 < m < 2000 GeV)";
    m_genxspb = 0.000038380;
    m_calxspb = 0.000038380;
    m_latex = "\\zrtautau[1750,2000]";
    m_run = 129553;
  } else if ( m_name == "ztautau_2000_2250" ) {
    m_title = "Z #rightarrow #tau #tau (2000 < m < 2250 GeV)";
    m_genxspb = 0.000013890;
    m_calxspb = 0.000013890;
    m_latex = "\\zrtautau[2000,2250]";
    m_run = 129554;
  } else if ( m_name == "ztautau_2250_2500" ) {
    m_title = "Z #rightarrow #tau #tau (2250 < m < 2500 GeV)";
    m_genxspb = 0.0000052260;
    m_calxspb = 0.0000052260;
    m_latex = "\\zrtautau[2250,2500]";
    m_run = 129555;
  } else if ( m_name == "ztautau_2500_2750" ) {
    m_title = "Z #rightarrow #tau #tau (2500 < m < 2750 GeV)";
    m_genxspb = 0.0000020170;
    m_calxspb = 0.0000020170;
    m_latex = "\\zrtautau[2500,2750]";
    m_run = 129556;
  } else if ( m_name == "ztautau_2750_3000" ) {
    m_title = "Z #rightarrow #tau #tau (2750 < m < 3000 GeV)";
    m_genxspb = 0.00000078910;
    m_calxspb = 0.00000078910;
    m_latex = "\\zrtautau[2750,3000]";
    m_run = 129557;
  } else if ( m_name == "ztautau_3000" ) {
    m_title = "Z #rightarrow #tau #tau (m > 3000 GeV)";
    m_genxspb = 0.00000050390;
    m_calxspb = 0.00000050390;
    m_latex = "\\rtautau[3000+]";
    m_run = 129558;
  // Z lepton pt filtered
  } else if ( m_name == "zee_pt100" ) {
    m_title = "Z #rightarrow ee (p_{T}^{e} > 100 GeV)";
    m_latex = "\\zee\\ ($\\pte > 100 \\gev$)";
    m_genxspb = 8.2464; // 41.230 * 0.20001
    m_calxspb = 8.2464;
    m_run = 129596;
  } else if ( m_name == "zee_pt150" ) {
    m_title = "Z #rightarrow ee (p_{T}^{e} > 150 GeV)";
    m_latex = "\\zee\\ ($\\pte > 150 \\gev$)";
    m_genxspb = 1.3095; //41.220 * 0.031769
    m_calxspb = 1.3095;
    m_run = 129626;
  } else if ( m_name == "zee_pt200" ) {
    m_title = "Z #rightarrow ee (p_{T}^{e} > 200 GeV)";
    m_latex = "\\zee\\ ($\\pte > 200 \\gev$)";
    m_genxspb = 0.41150; // 1.6510 * 0.24924
    m_calxspb = 0.41150;
    m_run = 129636;
  } else if ( m_name == "zee_pt300" ) {
    m_title = "Z #rightarrow ee (p_{T}^{e} > 300 GeV)";
    m_latex = "\\zee\\ ($\\pte > 300 \\gev$)";
    m_genxspb = 0.067495; // 1.6510 * 0.040881
    m_calxspb = 0.067495;
    m_run = 129646;
  } else if ( m_name == "zmumu_pt100" ) {
    m_title = "Z #rightarrow #mu #mu (p_{T}^{#mu} > 100 GeV)";
    m_latex = "\\zmumu\\ ($\\ptmu > 100 \\gev$)";
    m_genxspb = 8.5536; // 41.220 * 0.20751
    m_calxspb = 8.5536;
    m_run = 129597;
  } else if ( m_name == "zmumu_pt150" ) {
    m_title = "Z #rightarrow #mu #mu (p_{T}^{#mu} > 150 GeV)";
    m_latex = "\\zmumu\\ ($\\ptmu > 150 \\gev$)";
    m_genxspb = 1.3577; // 41.220 * 0.032939
    m_calxspb = 1.3577;
    m_run = 129627;
  } else if ( m_name == "zmumu_pt200" ) {
    m_title = "Z #rightarrow #mu #mu (p_{T}^{#mu} > 200 GeV)";
    m_latex = "\\zmumu\\ ($\\ptmu > 200 \\gev$)";
    m_genxspb = 0.42776; // 1.6510 * 0.25909
    m_calxspb = 0.42776;
    m_run = 129637;
  } else if ( m_name == "zmumu_pt300" ) {
    m_title = "Z #rightarrow #mu #mu (p_{T}^{#mu} > 300 GeV)";
    m_latex = "\\zmumu\\ ($\\ptmu > 300 \\gev$)";
    m_genxspb = 0.070742; // 1.6510 * 0.042848
    m_calxspb = 0.070742;
    m_run = 129647;
  // W/Z
  } else if ( m_name == "wzpt100_mu" ) {
     m_title = "W/Z #rightarrow #mu X (p_{T}^{#mu} > 100 GeV)";
     m_latex = "\\wzbg\\ ($\\ptmu > 100 \\gev$)";
     m_genxspb = 39.7;//35.7
     m_calxspb = 46.4;//41.7  // Using the K-factor for W->lnu (1.17)
     m_run = 119267;

  } else if ( m_name == "wzpt300_mu" ) {
    m_title = "W/Z #rightarrow #mu X (p_{T}^{#mu} > 300 GeV)";
    m_latex = "\\wzbg\\ ($\\ptmu > 300 \\gev$)";
    m_genxspb = 0.246;
    m_calxspb = 0.288; // Using the K-factor for W->lnu (1.17)
    m_run = 119607;
  } else if ( m_name == "wzpt100_e" ) {
    m_title = "W/Z #rightarrow e X (p_{T}^{e} > 100 GeV)";
    m_latex = "\\wzbg\\ ($\\ptmu > 100 \\gev$)";
    m_genxspb = 38.2;
    m_calxspb = 44.7; // Using the K-factor for W->lnu (1.17) 
    m_run = 119605;
    
  } else if ( m_name == "wzpt300_e" ) {
    m_title = "W/Z #rightarrow e X (p_{T}^{e} > 300 GeV)";
    m_latex = "\\wzbg\\ ($\\ptmu > 300 \\gev$)";
    m_genxspb = 0.236;
    m_calxspb = 0.276; // Using the K-factor for W-> lnu (1.17)
    m_run = 119606;

  // Diboson.
  } else if ( m_name == "wwl_herw" ) {
     m_title = "W W #rightarrow l X";
     m_latex = "\\wwdib";
     m_genxspb = 12.416; //32.501pb x 0.38203 (filt.eff)
     m_calxspb = 20.9005; //k-factor 1.6833 change cs to 21.7104?
     m_run = 105985;
  } else if ( m_name == "wwl" ) {
     m_title = "W W #rightarrow l X";
     m_latex = "\\wwdib";
     m_genxspb = 5.499; 
     m_calxspb = 5.829; //k-factor 1.06 //5.884; //k-factor 1.07
     m_run = 126892;
  } else if ( m_name == "wzl_herw" ) {
     m_title = "W Z #rightarrow l X";
     m_latex = "\\wzdib";
     m_genxspb = 3.6666; //12.009pb x 0.30532 (filt.eff)
     m_calxspb = 6.9706; //k-factor 1.9011 change cs to 6.5576?
     m_run = 105987;
  } else if ( m_name == "wz_lllnu" ) {
     m_title = "W Z #rightarrow lll nu";
     m_latex = "\\wzdib";
     m_genxspb = 9.7508; 
     m_calxspb = 10.238; //k-factor 1.05 //10.336; //k-factor 1.06
     m_run = 126893;
  } else if ( m_name == "wz_lnununu" ) {
     m_title = "W Z #rightarrow l nununu";
     m_latex = "\\wzdib";
     m_genxspb = 1.4047;
     m_calxspb = 1.4749; //k-factor 1.05
     m_run = 179975;
  } else if ( m_name == "zzl_herw" ) {
     m_title = "Z Z #rightarrow l X";
     m_latex = "\\zzdib";
     m_genxspb = 0.99244; //4.6915pb x 0.21154 (filt.eff)
     m_calxspb = 1.5379; //k-factor 1.5496 change cs to 1.5566?
     m_run = 105986;
  } else if ( m_name == "zz_llll" ) {
     m_title = "Z Z #rightarrow ll ll";
     m_latex = "\\zzdib";
     m_genxspb = 8.7403;
     m_calxspb = 8.7403; //k-factor 1.00 //9.7017; //k-factor 1.11
     m_run = 126894;
  } else if ( m_name == "zz_llnunu" ) {
    m_title = "Z Z #rightarrow ll nunu";
    m_latex = "\\zzdib";
    m_genxspb = 0.49601;
    m_calxspb = 0.52081; //k-factor 1.05 //0.56545; //k-factor 1.14
    m_run = 126895;
  } else if ( m_name == "vv_eeqq" ) {
    m_title = "V V #rightarrow ee qq";
    m_latex = "\\vveeqq";
    m_genxspb = 1.7020;
    m_calxspb = 1.7020;
    m_run = 157814;
  } else if ( m_name == "vv_mumuqq" ) {
    m_title = "V V #rightarrow #mu #mu qq";
    m_latex = "\\vvmumuqq";
    m_genxspb = 1.6870;
    m_calxspb = 1.6870;
    m_run = 157815;
  } else if ( m_name == "vv_tautauqq" ) {
    m_title = "V V #rightarrow #tau #tau qq";
    m_latex = "\\vvtautauqq";
    m_genxspb = 1.7023;
    m_calxspb = 1.7023;
    m_run = 157816;
  } else if ( m_name == "vv_enuqq" ) {
    m_title = "V V #rightarrow e #nu qq";
    m_latex = "\\vvenuqq";
    m_genxspb = 9.5566;
    m_calxspb = 9.5566;
    m_run = 157817;
  } else if ( m_name == "vv_munuqq" ) {
    m_title = "V V #rightarrow #mu #nu qq";
    m_latex = "\\vvmunuqq";
    m_genxspb = 9.5544;
    m_calxspb = 9.5544;
    m_run = 157818;
  } else if ( m_name == "vv_taunuqq" ) {
    m_title = "V V #rightarrow #tau #nu qq";
    m_latex = "\\vvtaunuqq";
    m_genxspb = 9.5568;
    m_calxspb = 9.5568;
    m_run = 157819;
  } else if (m_name == "wgmu" ) {
    m_title = "W #gamma #rightarrow #mu #nu #gamma";
    m_latex = "\\wgmudib";
    m_genxspb = 72.873;//7.2873E-02 nb 
    m_calxspb = 75.5;
    m_run = 126014;
  } else if ( m_name == "wge" ) {
    m_title = "W #gamma #rightarrow e #nu #gamma";
    m_latex = "\\wgedib";
    m_genxspb = 72.870; // 7.2870E-02 nb
    m_calxspb = 75.5;
    m_run = 126013;
  } else if ( m_name == "wgl_np0" ) {
    m_title = "W #gamma #rightarrow l #nu #gamma + 0 jets";
    m_latex = "\\wgldib_np0";
    m_genxspb = 72.118; //229.88pb * 0.31372 filt. eff.
    m_calxspb = 82.936; //k-factor 1.15
    m_run = 146436;
  } else if ( m_name == "wgl_np1" ) {
    m_title = "W #gamma #rightarrow l #nu #gamma + 1 jets";
    m_latex = "\\wgldib_np1";
    m_genxspb = 26.706; //59.518pb * 0.44871 filt. eff.
    m_calxspb = 30.712; //k-factor 1.15
    m_run = 146437;
  } else if ( m_name == "wgl_np2" ) {
    m_title = "W #gamma #rightarrow l #nu #gamma + 2 jets";
    m_latex = "\\wgldib_np2";
    m_genxspb = 11.649; //21.390pb * 0.54461 filt. eff.
    m_calxspb = 13.396; //k-factor 1.15
    m_run = 146438;
  } else if ( m_name == "wgl_np3" ) {
    m_title = "W #gamma #rightarrow l #nu #gamma + 3 jets";
    m_latex = "\\wgldib_np3";
    m_genxspb = 4.4839; //7.1203pb * 0.62974 filt. eff.
    m_calxspb = 5.156; //k-factor 1.15
    m_run = 146439;

  // ttbar
  } else if ( m_name == "ttbar_old" ) {
     m_title = "t tbar #rightarrow l X";
     m_latex = "\\top";
     m_genxspb = 112.936;//208.13pb x 0.54262 filt. eff.
     m_calxspb = 129.312;//k-factor of 1.129
     m_run = 105200;
  } else if ( m_name == "ttbar" ) {
     m_title = "t tbar #rightarrow l X";
     m_latex = "\\top";
     m_genxspb = 114.464;//210.83pb x 0.54292 filt. eff.
     m_calxspb = 137.32;//k-factor of 1.1997  //129.230;//k-factor of 1.129
     m_run = 105861;
  // Sigle Top
  } else if ( m_name == "st_wt" ) {
    m_title = "Wt inclusive";
    m_latex = "\\Wt";
    m_genxspb = 20.665;//
    m_calxspb = 22.380;//k-factor of 1.083
    m_run = 108346;
  } else if ( m_name == "st_tchan_enu" ) {
    m_title = "st_tchan_enu";
    m_latex = "\\t (t channel) #rightarrow enu";
    m_genxspb = 8.588;//
    m_calxspb =9.481;//k-factor of 1.104
    m_run = 117360;
  } else if ( m_name == "st_tchan_munu" ) {
    m_title = "st_tchan_munu";
    m_latex = "\\t (t channel) #rightarrow munu";
    m_genxspb = 8.589;
    m_calxspb = 9.482;//k-factor 1.104
    m_run = 117361;
  } else if ( m_name == "st_tchan_taunu" ) {
    m_title = "st_tchan_taunu";
    m_latex = "\\t (t channel #rightarrow taunu";
    m_genxspb = 8.581;
    m_calxspb = 9.482;//k-factor 1.105
    m_run = 117362;
  } else if ( m_name == "st_schan_enu" ) {
    m_title = "st_schan_enu";
    m_latex = "\\t (s channel #rightarrow enu";
    m_genxspb = 0.56394;
    m_calxspb = 0.606;//k-factor 1.074
    m_run = 108343;
  } else if ( m_name == "st_schan_munu" ) {
    m_title = "st_schan_munu";
    m_latex = "\\t (s channel #rightarrow munu";
    m_genxspb = 0.56431;//
    m_calxspb = 0.606;//k-factor 1.074
    m_run = 108344;
  } else if ( m_name == "st_schan_taunu" ) {
    m_title = "st_schan_taunu";
    m_latex = "\\t (s channel #rightarrow taunu";
    m_genxspb = 0.56436;//
    m_calxspb = 0.606;//k-factor 1.074
    m_run = 108345;

  // Wprime
  } else if ( m_name == "wpl150" ) {
     m_title = "W'(150) #rightarrow l #nu";
     m_latex = "\\wplm[150]";
     m_genxspb = 25.668;//25.254;//1299.7;
     m_calxspb = 25.668;//25.254;//1299.7;
     m_run = 182599; //158781; //115108;
  } else if ( m_name == "wpl150_old" ) {
     m_title = "W'(150) #rightarrow l #nu";
     m_latex = "\\wplm[150]";
     m_genxspb = 25.254;//1299.7;
     m_calxspb = 25.254;//1299.7;
     m_run = 158781; //115108;
  } else if ( m_name == "wpl200" ) {
     m_title = "W'(200) #rightarrow l #nu";
     m_latex = "\\wplm[200]";
     m_genxspb = 25.668;//25.254;//494.65;
     m_calxspb = 25.668;//25.254;//494.65;
     m_run = 182599; //158781; //115109;
  } else if ( m_name == "wpl200_old" ) {
     m_title = "W'(200) #rightarrow l #nu";
     m_latex = "\\wplm[200]";
     m_genxspb = 25.254;//494.65;
     m_calxspb = 25.254;//494.65;
     m_run = 158781; //115109;
  } else if ( m_name == "wpl300" ) {
     m_title = "W'(300) #rightarrow l #nu";
     m_latex = "\\wplm[300]";
     m_genxspb = 25.668;//25.254;//107.53; //109.35;
     m_calxspb = 157.74;//130.50; //109.35;
     m_run = 182599; //158781; //115110;
  } else if ( m_name == "wpl300_old" ) {
     m_title = "W'(300) #rightarrow l #nu";
     m_latex = "\\wplm[300]";
     m_genxspb = 25.254;//107.53; //109.35;
     m_calxspb = 157.74;//130.50; //109.35;
     m_run = 158781; //115110;
  } else if ( m_name == "wpl400" ) {
     m_title = "W'(400) #rightarrow l #nu";
     m_latex = "\\wplm[400]";
     m_genxspb = 25.668;//25.254;
     m_calxspb = 52.06; //53.78;
     m_run = 182599; //158781; //115111;
  } else if ( m_name == "wpl400_old" ) {
     m_title = "W'(400) #rightarrow l #nu";
     m_latex = "\\wplm[400]";
     m_genxspb = 25.254;
     m_calxspb = 52.06; //53.78;
     m_run = 158781; //115111;
  } else if ( m_name == "wpl500" ) {
     m_title = "W'(500) #rightarrow l #nu";
     m_latex = "\\wplm[500]";
     m_genxspb = 25.668;//25.254;
     m_calxspb = 22.02; //22.89;
     m_run = 182599; //158781; //115100;
  } else if ( m_name == "wpl500_old" ) {
     m_title = "W'(500) #rightarrow l #nu";
     m_latex = "\\wplm[500]";
     m_genxspb = 25.254;
     m_calxspb = 22.02; //22.89;
     m_run = 158781; //115100;
  } else if ( m_name == "wpl500_fix" ) {
    m_title = "W'(500) #rightarrow l #nu";
    m_latex = "\\wplm[500]";
    m_genxspb = 50.152; //15.55;//45.875 from ami / 2.9502
    m_calxspb = 22.02; //22.89; 
    m_run = 158760;
  } else if ( m_name == "wpl600" ) {
     m_title = "W'(600) #rightarrow l #nu";
     m_latex = "\\wplm[600]";
     m_genxspb = 25.668;//25.254;
     m_calxspb = 10.78688;//10.83; //11.17;
     m_run = 182599; //158781; //119657;
  } else if ( m_name == "wpl600_old" ) {
     m_title = "W'(600) #rightarrow l #nu";
     m_latex = "\\wplm[600]";
     m_genxspb = 25.254;
     m_calxspb = 10.78688;//10.83; //11.17;
     m_run = 158781; //119657;
  } else if ( m_name == "wpl750" ) {
     m_title = "W'(750) #rightarrow l #nu";
     m_latex = "\\wplm[750]";
     m_genxspb = 25.668;//25.254;//2.843;//2.931;
     m_calxspb = 4.30934; //4.429; //4.460;
     m_run = 182599; //158781; //115373;
  } else if ( m_name == "wpl750_old" ) {
     m_title = "W'(750) #rightarrow l #nu";
     m_latex = "\\wplm[750]";
     m_genxspb = 25.254;//2.843;//2.931;
     m_calxspb = 4.30934; //4.429; //4.460;
     m_run = 158781; //115373;
  } else if ( m_name == "wpl1000" ) {
     m_title = "W'(1000) #rightarrow l #nu";
     m_latex = "\\wplm[1000]";
     m_genxspb = 25.668;//25.254;//0.7659;//0.7792;
     m_calxspb = 1.204142;//1.257; //1.247;
     m_run = 182599; //158781; //115610;
  } else if ( m_name == "wpl1000_old" ) {
     m_title = "W'(1000) #rightarrow l #nu";
     m_latex = "\\wplm[1000]";
     m_genxspb = 25.254;//0.7659;//0.7792;
     m_calxspb = 1.204142;//1.257; //1.247;
     m_run = 158781; //115610;
  } else if ( m_name == "wpl1000_fix" ) {
     m_title = "W'(1000) #rightarrow l #nu";
     m_latex = "\\wplm[1000]";
     m_genxspb = 2.81;//0.7792;
     m_calxspb = 1.204142;//1.257; //1.247;
     m_run = 158761;
  } else if ( m_name == "wpl1250" ) {
     m_title = "W'(1250) #rightarrow l #nu";
     m_latex = "\\wplm[1250]";
     m_genxspb = 25.668;//25.254;//0.2451;//0.2555;
     m_calxspb = 0.403819; //0.4141; //0.4179;
     m_run = 182599; //158781; //105277;
  } else if ( m_name == "wpl1250_old" ) {
     m_title = "W'(1250) #rightarrow l #nu";
     m_latex = "\\wplm[1250]";
     m_genxspb = 25.254;//0.2451;//0.2555;
     m_calxspb = 0.403819; //0.4141; //0.4179;
     m_run = 158781; //105277;
  } else if ( m_name == "wpl1500" ) {
     m_title = "W'(1500) #rightarrow l #nu";
     m_latex = "\\wplm[1500]";
     m_genxspb = 25.668;//25.254;//0.08777;//0.09257;
     m_calxspb = 0.1502263; //0.1676; //0.1554;
     m_run = 182599; //158781; //115626;
  } else if ( m_name == "wpl1500_old" ) {
     m_title = "W'(1500) #rightarrow l #nu";
     m_latex = "\\wplm[1500]";
     m_genxspb = 25.254;//0.08777;//0.09257;
     m_calxspb = 0.1502263; //0.1676; //0.1554;
     m_run = 158781; //115626;
  } else if ( m_name == "wpl1750" ) {
     m_title = "W'(1750) #rightarrow l #nu";
     m_latex = "\\wplm[1750]";
     m_genxspb = 25.668;//25.254;//0.0337;//0.03622;
     m_calxspb = 0.059713; //0.05727; //0.0617;
     m_run = 182599; //158781; //105278;
  } else if ( m_name == "wpl1750_old" ) {
     m_title = "W'(1750) #rightarrow l #nu";
     m_latex = "\\wplm[1750]";
     m_genxspb = 25.254;//0.0337;//0.03622;
     m_calxspb = 0.059713; //0.05727; //0.0617;
     m_run = 158781; //105278;
  } else if ( m_name == "wpl2000" ) {
     m_title = "W'(2000) #rightarrow l #nu";
     m_latex = "\\wplm[2000]";
     m_genxspb = 25.668;//25.254;//0.01326;//0.01503;
     m_calxspb = 0.0248399; //0.02797; //0.0257;
     m_run = 182599; //158781; //105264;
  } else if ( m_name == "wpl2000_old" ) {
     m_title = "W'(2000) #rightarrow l #nu";
     m_latex = "\\wplm[2000]";
     m_genxspb = 25.254;//0.01326;//0.01503;
     m_calxspb = 0.0248399; //0.02797; //0.0257;
     m_run = 158781; //105264;
  } else if ( m_name == "wpl2250" ) {
     m_title = "W'(2250) #rightarrow l #nu";
     m_latex = "\\wplm[2250]";
     m_genxspb = 25.668;//25.254;//0.00605;//0.006618;
     m_calxspb = 0.0108098; //0.01082; //0.0112;
     m_run = 182599; //158781; //119658;
  } else if ( m_name == "wpl2250_old" ) {
     m_title = "W'(2250) #rightarrow l #nu";
     m_latex = "\\wplm[2250]";
     m_genxspb = 25.254;//0.00605;//0.006618;
     m_calxspb = 0.0108098; //0.01082; //0.0112;
     m_run = 158781; //119658;
  } else if ( m_name == "wpl2500" ) {
     m_title = "W'(2500) #rightarrow l #nu";
     m_latex = "\\wplm[2500]";
     m_genxspb = 25.668;//25.254;//0.00289;//0.00314;
     m_calxspb = 0.004933; //0.004905; //0.0051;
     m_run = 182599; //158781; //105265;
  } else if ( m_name == "wpl2500_old" ) {
     m_title = "W'(2500) #rightarrow l #nu";
     m_latex = "\\wplm[2500]";
     m_genxspb = 25.254;//0.00289;//0.00314;
     m_calxspb = 0.004933; //0.004905; //0.0051;
     m_run = 158781; //105265;
  } else if ( m_name == "wpl2750" ) {
     m_title = "W'(2750) #rightarrow l #nu";
     m_latex = "\\wplm[2750]";
     m_genxspb = 25.668;//25.254;//0.001559;//0.00154;
     m_calxspb = 0.0024027; //0.002436; //0.0025;
     m_run = 182599; //158781; //119659;
  } else if ( m_name == "wpl2750_old" ) {
     m_title = "W'(2750) #rightarrow l #nu";
     m_latex = "\\wplm[2750]";
     m_genxspb = 25.254;//0.001559;//0.00154;
     m_calxspb = 0.0024027; //0.002436; //0.0025;
     m_run = 158781; //119659;
  } else if ( m_name == "wpl3000" ) {
     m_title = "W'(3000) #rightarrow l #nu";
     m_latex = "\\wplm[3000]";
     m_genxspb = 25.668;//25.254;//0.000872;//0.000902;
     m_calxspb = 0.0012662; //0.001297; //0.0013;
     m_run = 182599; //158781; //105663;
  } else if ( m_name == "wpl3000_old" ) {
     m_title = "W'(3000) #rightarrow l #nu";
     m_latex = "\\wplm[3000]";
     m_genxspb = 25.254;//0.000872;//0.000902;
     m_calxspb = 0.0012662; //0.001297; //0.0013;
     m_run = 158781; //105663;
  } else if ( m_name == "wpl3000_fix" ) {
     m_title = "W'(3000) #rightarrow l #nu";
     m_latex = "\\wplm[3000]";
     m_genxspb = 0.0038665;
     m_calxspb = 0.0012662; //0.001297; //0.0013;
     m_run = 158762;
  } else if ( m_name == "wpl3250" ) {
     m_title = "W'(3250) #rightarrow l #nu";
     m_latex = "\\wplm[3250]";
     m_genxspb = 25.668;//25.254;//0.000872;//0.000902;
     m_calxspb = 0.0007288; //0.0007316; //0.0008;
     m_run = 182599; //158781; //105663;
  } else if ( m_name == "wpl3250_old" ) {
     m_title = "W'(3250) #rightarrow l #nu";
     m_latex = "\\wplm[3250]";
     m_genxspb = 25.254;//0.000872;//0.000902; 
     m_calxspb = 0.0007288; //0.0007316; //0.0008;
     m_run = 158781; //105663;
  } else if ( m_name == "wpl3500" ) {
     m_title = "W'(3500) #rightarrow l #nu";
     m_latex = "\\wplm[3500]";
     m_genxspb = 25.668;//25.254;//0.000388;//0.000393;
     m_calxspb = 0.0004566; //0.0004578; //0.0005;
     m_run = 182599; //158781; //105664;
  } else if ( m_name == "wpl3500_old" ) {
     m_title = "W'(3500) #rightarrow l #nu";
     m_latex = "\\wplm[3500]";
     m_genxspb = 25.254;//0.000388;//0.000393;
     m_calxspb = 0.0004566; //0.0004578; //0.0005;
     m_run = 158781; //105664;
  } else if ( m_name == "wpl3750" ) {
     m_title = "W'(3750) #rightarrow l #nu";
     m_latex = "\\wplm[3750]";
     m_genxspb = 25.668;//25.254;//0.000388;//0.000902;
     m_calxspb = 0.0003061; //0.000312; //0.0003;//0.000221;//0.000902;
     m_run = 182599; //158781; //105663;
  } else if ( m_name == "wpl3750_old" ) {
     m_title = "W'(3750) #rightarrow l #nu";
     m_latex = "\\wplm[3750]";
     m_genxspb = 25.254;//0.000388;//0.000902; 
     m_calxspb = 0.0003061; //0.000312; //0.0003;//0.000221;//0.000902; 
     m_run = 158781; //105663;
  } else if ( m_name == "wpl4000" ) {
     m_title = "W'(4000) #rightarrow l #nu";
     m_latex = "\\wplm[4000]";
     m_genxspb = 25.668;//25.254;//0.000388;//0.000902;
     m_calxspb = 0.0002163; //0.0002171; //0.0002;
     m_run = 182599; //158781;
  } else if ( m_name == "wpl4000_old" ) {
     m_title = "W'(4000) #rightarrow l #nu";
     m_latex = "\\wplm[4000]";
     m_genxspb = 25.254;//0.000388;//0.000902; 
     m_calxspb = 0.0002163; //0.0002171; //0.0002;
     m_run = 158781;

  //W star samples
 } else if ( m_name == "wse150" ) {
     m_title = "W*(150) #rightarrow e #nu";
     m_latex = "\\wsem[150]";
     m_genxspb = 8.455;
     m_calxspb = 8.455;
     m_run = 158763;
 } else if ( m_name == "wse200" ) {
     m_title = "W*(200) #rightarrow e #nu";
     m_latex = "\\wsem[200]";
     m_genxspb = 8.455;
     m_calxspb = 8.455;
     m_run = 158763;
 } else if ( m_name == "wse300" ) {
     m_title = "W*(300) #rightarrow e #nu";
     m_latex = "\\wsem[300]";
     m_genxspb = 8.455;
     m_calxspb = 61.441;
     m_run = 158763;
 } else if ( m_name == "wse400" ) {
     m_title = "W*(400) #rightarrow e #nu";
     m_latex = "\\wsem[400]";
     m_genxspb = 8.455;//32.834;//29.59;
     m_calxspb = 37.907;//32.834;
     m_run = 158763;
  } else if ( m_name == "wse500" ) {
     m_title = "W*(500) #rightarrow e #nu";
     m_latex = "\\wsem[500]";
     m_genxspb = 8.455;//13.612;//12.60;//45.875 from ami / 2.9502
     m_calxspb = 16.451;//13.612;//12.60;
     m_run = 158763;
  } else if ( m_name == "wse500_fix" ) {
     m_title = "W*(500) #rightarrow e #nu";
     m_latex = "\\wsem[500]";
     m_genxspb = 16.22;//13.612;
     m_calxspb = 16.451;//13.612;
     m_run = 158751;
  } else if ( m_name == "wse600" ) {
     m_title = "W*(600) #rightarrow e #nu";
     m_latex = "\\wsem[600]";
     m_genxspb = 8.455;
     m_calxspb = 8.0757;
     m_run = 158763;
  } else if ( m_name == "wse750" ) {
     m_title = "W*(750) #rightarrow e #nu";
     m_latex = "\\wsem[750]";
     m_genxspb = 8.455;//2.4309;//2.34;
     m_calxspb = 3.2428;//2.4309;//2.34;
     m_run = 158763;
  } else if ( m_name == "wse1000" ) {
     m_title = "W*(1000) #rightarrow e #nu";
     m_latex = "\\wsem[1000]";
     m_genxspb = 8.455;//0.61184;//0.6102;
     m_calxspb = 0.90386;//0.61184;//0.6102;
     m_run = 158763;
  } else if ( m_name == "wse1000_fix" ) {
     m_title = "W*(1000) #rightarrow e #nu";
     m_latex = "\\wsem[1000]";
     m_genxspb = 0.8836;//0.61184;
     m_calxspb = 0.90386;//0.61184;
     m_run = 158752;
  } else if ( m_name == "wse1250" ) {
     m_title = "W*(1250) #rightarrow e #nu";
     m_latex = "\\wsem[1250]";
     m_genxspb = 8.455;//0.18194;//0.1883;
     m_calxspb = 0.29909;//0.18194;//0.1883;
     m_run = 158763;
  } else if ( m_name == "wse1500" ) {
     m_title = "W*(1500) #rightarrow e #nu";
     m_latex = "\\wsem[1500]";
     m_genxspb = 8.455;//0.059298;//0.06360;
     m_calxspb = 0.10885;//0.059298;//0.06360;
     m_run = 158763;
  } else if ( m_name == "wse1750" ) {
     m_title = "W*(1750) #rightarrow e #nu";
     m_latex = "\\wsem[1750]";
     m_genxspb = 8.455;//0.020278;//0.02256;
     m_calxspb = 0.041734;//0.020278;//0.02256;
     m_run = 158763;
  } else if ( m_name == "wse2000" ) {
     m_title = "W*(2000) #rightarrow e #nu";
     m_latex = "\\wsem[2000]";
     m_genxspb = 8.455;//0.0071366;//0.00819;
     m_calxspb = 0.016505;//0.0071366;//0.00819;
     m_run = 158763;
  } else if ( m_name == "wse2250" ) {
     m_title = "W*(2250) #rightarrow e #nu";
     m_latex = "\\wsem[2250]";
     m_genxspb = 8.455;//0.0025477;//0.002994;
     m_calxspb = 0.0066168;//0.0025477;//0.002994;
     m_run = 158763;
  } else if ( m_name == "wse2500" ) {
     m_title = "W*(2500) #rightarrow e #nu";
     m_latex = "\\wsem[2500]";
     m_genxspb = 8.455;//0.00091657;//0.001088;
     m_calxspb = 0.0026713;//0.00091657;//0.001088;
     m_run = 158763;
  } else if ( m_name == "wse2750" ) {
     m_title = "W*(2750) #rightarrow e #nu";
     m_latex = "\\wsem[2750]";
     m_genxspb = 8.455;//0.00033089;//0.000391;
     m_calxspb = 0.0010801;//0.00033089;//0.000391;
     m_run = 158763;
  } else if ( m_name == "wse3000" ) {
     m_title = "W*(3000) #rightarrow e #nu";
     m_latex = "\\wsem[3000]";
     m_genxspb = 8.455;//0.00012051;//0.000138;
     m_calxspb = 0.00043621;//0.00012051;//0.000138;
     m_run = 158763;
  } else if ( m_name == "wse3000_fix" ) {
     m_title = "W*(3000) #rightarrow e #nu";
     m_latex = "\\wsem[3000]";
     m_genxspb = 0.0004896;//0.00012051;
     m_calxspb = 0.00043621;//0.00012051;
     m_run = 158753;
  } else if ( m_name == "wse3250" ) {
     m_title = "W*(3250) #rightarrow e #nu";
     m_latex = "\\wsem[3250]";
     m_genxspb = 8.455;
     m_calxspb = 0.00017587;
     m_run = 158763;
  } else if ( m_name == "wse3500" ) {
     m_title = "W*(3500) #rightarrow e #nu";
     m_latex = "\\wsem[3500]";
     m_genxspb = 8.455;
     m_calxspb = 0.000071277;
     m_run = 158763;
  } else if ( m_name == "wse3750" ) {
     m_title = "W*(3750) #rightarrow e #nu";
     m_latex = "\\wsem[3750]";
     m_genxspb = 8.455;
     m_calxspb = 0.000029122;
     m_run = 158763;
  } else if ( m_name == "wse4000" ) {
     m_title = "W*(4000) #rightarrow e #nu";
     m_latex = "\\wsem[4000]";
     m_genxspb = 8.455;
     m_calxspb = 0.000012143;
     m_run = 158763;

 } else if ( m_name == "wsm150" ) {
     m_title = "W*(150) #rightarrow mu #nu";
     m_latex = "\\wsmm[150]";
     m_genxspb = 8.449;
     m_calxspb = 8.449;
     m_run = 158764;
 } else if ( m_name == "wsm200" ) {
     m_title = "W*(200) #rightarrow mu #nu";
     m_latex = "\\wsmm[200]";
     m_genxspb = 8.449;
     m_calxspb = 8.449;
     m_run = 158764;
 } else if ( m_name == "wsm300" ) {
     m_title = "W*(300) #rightarrow mu #nu";
     m_latex = "\\wsmm[300]";
     m_genxspb = 8.449;
     m_calxspb = 61.441;
     m_run = 158764;
 } else if ( m_name == "wsm400" ) {
     m_title = "W*(400) #rightarrow mu #nu";
     m_latex = "\\wsmm[400]";
     m_genxspb = 8.449;//32.834;//29.59;
     m_calxspb = 37.907;//32.834;//29.59;
     m_run = 158764;
  } else if ( m_name == "wsm500" ) {
     m_title = "W*(500) #rightarrow mu #nu";
     m_latex = "\\wsmm[500]";
     m_genxspb = 8.449;//13.612;//12.60;//45.875 from ami / 2.9502
     m_calxspb = 16.451;//13.612;//12.60;
     m_run = 158764;
  } else if ( m_name == "wsm500_fix" ) {
     m_title = "W*(500) #rightarrow mu #nu";
     m_latex = "\\wsmm[500]";
     m_genxspb = 16.26;//13.612;
     m_calxspb = 16.451;//13.612;
     m_run = 158754;
  } else if ( m_name == "wsm600" ) {
     m_title = "W*(600) #rightarrow mu #nu";
     m_latex = "\\wsmm[600]";
     m_genxspb = 8.449;
     m_calxspb = 8.0757;
     m_run = 158764;
  } else if ( m_name == "wsm750" ) {
     m_title = "W*(750) #rightarrow mu #nu";
     m_latex = "\\wsmm[750]";
     m_genxspb = 8.449;//2.4309;//2.34;
     m_calxspb = 3.2428;//2.4309;//2.34;
     m_run = 158764;
  } else if ( m_name == "wsm1000" ) {
     m_title = "W*(1000) #rightarrow mu #nu";
     m_latex = "\\wsmm[1000]";
     m_genxspb = 8.449;//0.61184;//0.6102;
     m_calxspb = 0.90386;//0.61184;//0.6102;
     m_run = 158764;
  } else if ( m_name == "wsm1000_fix" ) {
     m_title = "W*(1000) #rightarrow mu #nu";
     m_latex = "\\wsmm[1000]";
     m_genxspb = 0.8837;//0.61184;
     m_calxspb = 0.90386;//0.61184;
     m_run = 158755;
  } else if ( m_name == "wsm1250" ) {
     m_title = "W*(1250) #rightarrow mu #nu";
     m_latex = "\\wsmm[1250]";
     m_genxspb = 8.449;//0.18194;//0.1883;
     m_calxspb = 0.29909;//0.18194;//0.1883;
     m_run = 158764;
  } else if ( m_name == "wsm1500" ) {
     m_title = "W*(1500) #rightarrow mu #nu";
     m_latex = "\\wsmm[1500]";
     m_genxspb = 8.449;//0.059298;//0.06360;
     m_calxspb = 0.10885;//0.059298;//0.06360;
     m_run = 158764;
  } else if ( m_name == "wsm1750" ) {
     m_title = "W*(1750) #rightarrow mu #nu";
     m_latex = "\\wsmm[1750]";
     m_genxspb = 8.449;//0.020278;//0.02256;
     m_calxspb = 0.041734;//0.020278;//0.02256;
     m_run = 158764;
  } else if ( m_name == "wsm2000" ) {
     m_title = "W*(2000) #rightarrow mu #nu";
     m_latex = "\\wsmm[2000]";
     m_genxspb = 8.449;//0.0071366;//0.00819;
     m_calxspb = 0.016505;//0.0071366;//0.00819;
     m_run = 158764;
  } else if ( m_name == "wsm2250" ) {
     m_title = "W*(2250) #rightarrow mu #nu";
     m_latex = "\\wsmm[2250]";
     m_genxspb = 8.449;//0.0025477;//0.002994;
     m_calxspb = 0.0066168;//0.0025477;//0.002994;
     m_run = 158764;
  } else if ( m_name == "wsm2500" ) {
     m_title = "W*(2500) #rightarrow mu #nu";
     m_latex = "\\wsmm[2500]";
     m_genxspb = 8.449;//0.00091657;//0.001088;
     m_calxspb = 0.0026713;//0.00091657;//0.001088;
     m_run = 158764;
  } else if ( m_name == "wsm2750" ) {
     m_title = "W*(2750) #rightarrow mu #nu";
     m_latex = "\\wsmm[2750]";
     m_genxspb = 8.449;//0.00033089;//0.000391;
     m_calxspb = 0.0010801;//0.00033089;//0.000391;
     m_run = 158764;
  } else if ( m_name == "wsm3000" ) {
     m_title = "W*(3000) #rightarrow mu #nu";
     m_latex = "\\wsmm[3000]";
     m_genxspb = 8.449;//0.00012051;//0.000138;
     m_calxspb = 0.00043621;//0.00012051;//0.000138;
     m_run = 158764;
  } else if ( m_name == "wsm3000_fix" ) {
     m_title = "W*(3000) #rightarrow mu #nu";
     m_latex = "\\wsmm[3000]";
     m_genxspb = 0.0004904;//0.00012051;
     m_calxspb = 0.00043621;//0.00012051;
     m_run = 158756;
  } else if ( m_name == "wsm3250" ) {
     m_title = "W*(3250) #rightarrow mu #nu";
     m_latex = "\\wsmm[3250]";
     m_genxspb = 8.449;
     m_calxspb = 0.00017587;
     m_run = 158764;
  } else if ( m_name == "wsm3500" ) {
     m_title = "W*(3500) #rightarrow mu #nu";
     m_latex = "\\wsmm[3500]";
     m_genxspb = 8.449;
     m_calxspb = 0.000071277;
     m_run = 158764;
  } else if ( m_name == "wsm3750" ) {
     m_title = "W*(3750) #rightarrow mu #nu";
     m_latex = "\\wsmm[3750]";
     m_genxspb = 8.449;
     m_calxspb = 0.000029122;
     m_run = 158764;
  } else if ( m_name == "wsm4000" ) {
     m_title = "W*(4000) #rightarrow mu #nu";
     m_latex = "\\wsmm[4000]";
     m_genxspb = 8.449;
     m_calxspb = 0.000012143;
     m_run = 158764;

  //mono-W samples
  } else if ( m_name == "d1_dm1" ) {
     m_title = "mono-W D1 (1)";
     m_latex = "\\d1[1]";
     m_genxspb = 4.3860E-10;
     m_calxspb = 4.3860E-10;
     m_run = 158820;
  } else if ( m_name == "d1_dm100" ) {
     m_title = "mono-W D1 (100)";
     m_latex = "\\d1[100]";
     m_genxspb = 3.3190E-10;
     m_calxspb = 3.3190E-10;
     m_run = 158821;
  } else if ( m_name == "d1_dm200" ) {
     m_title = "mono-W D1 (200)";
     m_latex = "\\d1[200]";
     m_genxspb = 2.0120E-10;
     m_calxspb = 2.0120E-10;
     m_run = 158822;
  } else if ( m_name == "d1_dm400" ) {
     m_title = "mono-W D1 (400)";
     m_latex = "\\d1[400]";
     m_genxspb = 6.4560E-11;
     m_calxspb = 6.4560E-11;
     m_run = 158823;
  } else if ( m_name == "d1_dm1000" ) {
     m_title = "mono-W D1 (1000)";
     m_latex = "\\d1[1000]";
     m_genxspb = 1.6030E-12;
     m_calxspb = 1.6030E-12;
     m_run = 158824;
  } else if ( m_name == "d1_dm1300" ) {
     m_title = "mono-W D1 (1300)";
     m_latex = "\\d1[1300]";
     m_genxspb = 2.1260E-13;
     m_calxspb = 2.1260E-13;
     m_run = 158825;
  } else if ( m_name == "d9_dm1" ) {
     m_title = "mono-W D9 (1)";
     m_latex = "\\d9[1]";
     m_genxspb = 9.6630E-02;
     m_calxspb = 9.6630E-02;
     m_run = 158826;
  } else if ( m_name == "d9_dm100") {
     m_title = "mono-W D9 (100)";
     m_latex = "\\d9[100]";
     m_genxspb = 8.6970E-02;
     m_calxspb = 8.6970E-02;
     m_run = 158827;
  } else if ( m_name == "d9_dm200") {
     m_title = "mono-W D9 (200)";
     m_latex = "\\d9[200]";
     m_genxspb = 6.9450E-02;
     m_calxspb = 6.9450E-02;
     m_run = 158828;
  } else if ( m_name == "d9_dm400") {
     m_title = "mono-W D9 (400)";
     m_latex = "\\d9[400]";
     m_genxspb = 3.6500E-02;
     m_calxspb = 3.6500E-02;
     m_run = 158829;
  } else if ( m_name == "d9_dm1000") {
     m_title = "mono-W D9 (1000)";
     m_latex = "\\d9[1000]";
     m_genxspb = 2.2700E-03;
     m_calxspb = 2.2700E-03;
     m_run = 158830;
  } else if ( m_name == "d9_dm1300") {
     m_title = "mono-W D9 (1300)";
     m_latex = "\\d9[1300]";
     m_genxspb = 4.1180E-04;
     m_calxspb = 4.1180E-04;
     m_run = 158831;
  } else if ( m_name == "d5_dm1") {
     m_title = "mono-W D5 (1)";
     m_latex = "\\d5[1]";
     m_genxspb = 7.2190E-03;
     m_calxspb = 7.2190E-03;
     m_run = 158832;
  } else if ( m_name == "d5_dm100") {
     m_title = "mono-W D5 (100)";
     m_latex = "\\d5[100]";
     m_genxspb = 7.0800E-03;
     m_calxspb = 7.0800E-03;
     m_run = 158833;
  } else if ( m_name == "d5_dm200") {
     m_title = "mono-W D5 (200)";
     m_latex = "\\d5[200]";
     m_genxspb = 5.8840E-03;
     m_calxspb = 5.8840E-03;
     m_run = 158834;
  } else if ( m_name == "d5_dm400") {
     m_title = "mono-W D5 (400)";
     m_latex = "\\d5[400]";
     m_genxspb = 3.2910E-03;
     m_calxspb = 3.2910E-03;
     m_run = 158835;
  } else if ( m_name == "d5_dm1000") {
     m_title = "mono-W D5 (1000)";
     m_latex = "\\d5[1000]";
     m_genxspb = 2.3690E-04;
     m_calxspb = 2.3690E-04;
     m_run = 158836;
  } else if ( m_name == "d5_dm1300") {
     m_title = "mono-W D5 (1300)";
     m_latex = "\\d5[1300]";
     m_genxspb = 4.5370E-05;
     m_calxspb = 4.5370E-05;
     m_run = 158837;
  } else if ( m_name == "d52_dm1") {
     m_title = "mono-W D52 (1)";
     m_latex = "\\d52[1]";
     m_genxspb = 6.0820E-02;
     m_calxspb = 6.0820E-02;
     m_run = 158838;
  } else if ( m_name == "d52_dm100") {
     m_title = "mono-W D52 (100)";
     m_latex = "\\d52[100]";
     m_genxspb = 5.7510E-02;
     m_calxspb = 5.7510E-02;
     m_run = 158839;
  } else if ( m_name == "d52_dm200") {
     m_title = "mono-W D52 (200)";
     m_latex = "\\d52[200]";
     m_genxspb = 4.8790E-02;
     m_calxspb = 4.8790E-02;
     m_run = 158840;
  } else if ( m_name == "d52_dm400") {
     m_title = "mono-W D52 (400)";
     m_latex = "\\d52[400]";
     m_genxspb = 2.7920E-02;
     m_calxspb = 2.7920E-02;
     m_run = 158841;
  } else if ( m_name == "d52_dm1000") {
     m_title = "mono-W D52 (1000)";
     m_latex = "\\d52[1000]";
     m_genxspb = 1.9170E-03;
     m_calxspb = 1.9170E-03;
     m_run = 158842;
  } else if ( m_name == "d52_dm1300") {
     m_title = "mono-W D52 (1300)";
     m_latex = "\\d52[1300]";
     m_genxspb = 3.5120E-04;
     m_calxspb = 3.5120E-04;
     m_run = 158843;


  // Data.
  } else if ( m_name == "data" ) {
     m_title = "Data";
     m_latex = "\\Data";
     m_genxspb = 0.;
     m_calxspb = 0.;
     m_run = 999;
  } else if ( m_name == "qcd" ) {
     m_title = "QCD";
     m_latex = "\\QCD";
     m_genxspb = 0.;
     m_calxspb = 0.;
     m_run = 9999;
  } else if ( m_name == "debug" ) {
     m_title = "Debug";
     m_latex = "\\Debug";
     m_genxspb = 0.0;
     m_calxspb = 0.0;
     m_run = 0;
  } 
  else {
    cerr << "WprimeSample::WprimeSample: Unknown sample: " << m_name << endl;
    m_name = "";
  }
}

std::ostream& operator<<(std::ostream& lhs, const WprimeSample& rhs) {
  lhs <<   "             Name: " << rhs.name();
  lhs << "\n            Alias: " << rhs.alias();
  lhs << "\n            Title: " << rhs.title();
  lhs << "\n   Generated xsec: " << rhs.genxspb() << " pb";
  lhs << "\n  Calculated xsec: " << rhs.calxspb() << " pb";
  lhs << "\n            Latex: " << rhs.latex();
  lhs << "\n              Run: " << rhs.run();
  return lhs;
}
