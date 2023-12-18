
#include <cinolib/meshes/meshes.h>
#include <cinolib/triangle_wrap.h>
using namespace std;
using namespace cinolib;
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <time.h>

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


void make_unit_circle(uint N, std::vector<double> & coords, std::vector<uint> & segs )
// generate a N-polygon approximating a unit circle centered at the origin
{
  coords.resize(0);
  segs.resize(0);
  for (uint i=0;i<N;i++) {
    coords.push_back(cos(i*2*M_PI/N)); 
    coords.push_back(sin(i*2*M_PI/N));
    segs.push_back(i);segs.push_back((i+1)%N);
  }
}

void make_domain(const std::vector<double> &coords_in, const std::vector<uint> &segs_in,
                  std::vector<double> &coords_out, std::vector<uint> &tris_out,
                  const double max_area, const double min_angle = 20) {
  std::string flags =
      "Qq" + std::to_string(min_angle) + "a" + std::to_string(max_area);
  triangle_wrap(coords_in, segs_in, {}, 0, flags.c_str(), coords_out, tris_out);
}

void make_disk(uint N,std::vector<double> &coords_out, std::vector<uint> &tris_out)
{
  std::vector<double> coords; 
  std::vector<uint> segs;
  make_unit_circle(N,coords,segs);
  make_domain(coords,segs,coords_out,tris_out,sqrt(3)*M_PI/(N*N),20);
}

void twist_disk(double delta, std::vector<vector<double>> &soup)
{
  uint n=soup.size();
  for (uint i=0;i<n;i++) 
    for (uint j=0;j<12;j+=2)
    {
      double x = soup[i][j], y = soup[i][j+1];
      double rho = sqrt(x*x+y*y);
      double theta = atan2(y,x)+delta*rho;
      soup[i][j] = rho * cos(theta);
      soup[i][j+1] = rho * sin(theta);
  }
}



void write_soup(const char *filename,const vector<vector<double>> &soup)
{
  ofstream f(filename);
  uint ntris = soup.size();
  f << ntris << endl;
  for (int i=0;i<ntris;i++ ) {
    for (uint j=0;j<12;j++) f << soup[i][j] << " ";
    f << endl;
  }
  f.close();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void make_soup(const std::vector<double> &coords, const std::vector<uint> &tris,
                std::vector<vector<double>> & soup)
{
  uint ntris = tris.size()/3;
  for (int i=0;i<ntris;i++) {
    soup[i][0]= coords[3*tris[3*i]];
    soup[i][1]= coords[3*tris[3*i]+1];
    soup[i][2]= coords[3*tris[3*i+1]];
    soup[i][3]= coords[3*tris[3*i+1]+1];
    soup[i][4]= coords[3*tris[3*i+2]];
    soup[i][5]= coords[3*tris[3*i+2]+1];
    soup[i][6]=(soup[i][0]+soup[i][2])/2.0;
    soup[i][7]=(soup[i][1]+soup[i][3])/2.0;
    soup[i][8]=(soup[i][2]+soup[i][4])/2.0;
    soup[i][9]=(soup[i][3]+soup[i][5])/2.0;
    soup[i][10]=(soup[i][4]+soup[i][0])/2.0;
    soup[i][11]=(soup[i][5]+soup[i][1])/2.0;
  }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                                 GUI utility

int main(int argc, char **argv) {
  std::vector<double> coords_out;
  std::vector<uint> tris_out;
  std::vector<uint> quads{};
  uint N = 10;  // default
  double delta = 10;  // rotation angle per step
  uint n_disks = 1;   // # steps-files
  if (argc > 1) N = std::stoi(argv[1]);
  if (argc > 2) n_disks = std::stoi(argv[2]);
  if (argc > 3) delta = std::stod(argv[3]);
  make_disk(N,coords_out,tris_out);
  std::vector<vector<double>> soup(tris_out.size()/3,vector<double>(12));
  make_soup(coords_out,tris_out,soup);
  std::string namebase = "../../output/disk_P2" + std::to_string(N);
  std::string name = namebase  + ".msh";
  write_soup(name.c_str(), soup);
  for (uint i=0;i<n_disks;i++) {
    twist_disk(delta*M_PI/180.0,soup);
    name = namebase + "_t" + std::to_string(static_cast<int>(delta*(i+1))) + ".msh";
    write_soup(name.c_str(), soup);
  }
}
