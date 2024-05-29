//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
# include "../headers/GnuplotC.h"
//////////////////////////////////////////////////////////
void GnuplotC::plot_scaling_relations(string prop)
{
  real_prec a=0;
}
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
void GnuplotC::plot_power_spectrum(vector<real_prec>&kv, vector<real_prec>&power)
{
  this->gp<<"set border linewidth 2.\n";
  this->gp<<"set size square \n";
  this->gp<<"set log x \n";
  this->gp<<"set xlabel 'k [h / Mpc]' font 'Times-Roman,14'\n";
  this->gp<<"set ylabel 'log P(k) [(Mpc / h)³]'  font 'Times-Roman,14'\n";
  vector<pair<real_prec, real_prec> > xy_pts;
  for(ULONG i=0; i<kv.size(); ++i)
    xy_pts.push_back(std::make_pair(kv[i], log10(power[i])));
  this->gp<<"plot"<<this->gp.file1d(xy_pts)<<" w l lw 3.3 lt 6 title 'Power'"<<endl;
}
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
void GnuplotC::plot_power_spectrum_redshift_space(vector<real_prec>&kv, vector<real_prec>&power0,vector<real_prec>&power2,vector<real_prec>&power4)
{
  this->gp<<"set border linewidth 2.\n";
  this->gp<<"set size square \n";
  this->gp<<"set log x \n";
  this->gp<<"set xlabel 'k [h / Mpc]' font 'Times-Roman,14'\n";
  this->gp<<"set ylabel 'log P(k) [(Mpc / h)³]'  font 'Times-Roman,14'\n";

  vector<pair<real_prec, real_prec> > xy_pts0;
  for(ULONG i=0; i<kv.size(); ++i)
    xy_pts0.push_back(std::make_pair(kv[i], log10(power0[i])));

  vector<pair<real_prec, real_prec> > xy_pts2;
  for(ULONG i=0; i<kv.size(); ++i)
    xy_pts2.push_back(std::make_pair(kv[i], log10(power2[i])));
  
  vector<pair<real_prec, real_prec> > xy_pts4;
  for(ULONG i=0; i<kv.size(); ++i)
    xy_pts4.push_back(std::make_pair(kv[i], log10(power4[i])));
  
  this->gp<<"plot"<<this->gp.file1d(xy_pts0)<<" w l lw 3. lt 1 title 'P0'"<<endl;
  this->gp<<"replot"<<this->gp.file1d(xy_pts2)<<" w l lw 3. lt 3 title 'P2'"<<endl;
  this->gp<<"replot"<<this->gp.file1d(xy_pts4)<<" w l lw 3. lt 4 title 'P4'"<<endl;
}
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
void GnuplotC::plot_power_spectrum_redshift_space(vector<real_prec>&kv, vector<real_prec>&power0,vector<real_prec>&power2)
{
  this->gp<<"set border linewidth 2.\n";
  this->gp<<"set size square \n";
  this->gp<<"set log x \n";
  this->gp<<"set xlabel 'k [h / Mpc]' font 'Times-Roman,14'\n";
  this->gp<<"set ylabel 'log P(k) [(Mpc / h)³]'  font 'Times-Roman,14'\n";

  vector<pair<real_prec, real_prec> > xy_pts0;
  for(ULONG i=0; i<kv.size(); ++i)
    xy_pts0.push_back(std::make_pair(kv[i], log10(power0[i])));

  vector<pair<real_prec, real_prec> > xy_pts2;
  for(ULONG i=0; i<kv.size(); ++i)
    xy_pts2.push_back(std::make_pair(kv[i], log10(power2[i])));


  this->gp<<"plot"<<this->gp.file1d(xy_pts0)<<" w l lw 3. lt 1 title 'l=0'"<<endl;
  this->gp<<"replot"<<this->gp.file1d(xy_pts2)<<" w l lw 3. lt 3 title 'l=2'"<<endl;
}
