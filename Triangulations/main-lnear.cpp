
#include <cinolib/gl/glcanvas.h>
#include <cinolib/io/write_OBJ.h>
#include <cinolib/meshes/meshes.h>
#include <cinolib/triangle_wrap.h>
using namespace std;
using namespace cinolib;
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <time.h>

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void regular_tri(std::vector<double> &coords_out, std::vector<uint> &tris_out,
                 int N) {
  double deltax = 1.0 / N;              // edge length
  double deltay = deltax * sqrt(3) / 2; // height of tris
  int Nh = (int)(1.0 / deltay);         // number of rows
  double shift = 0; // commented to stretch tris to y=1 (1.0-Nh*deltay)/2.0; //
                    // brings (0.5,0.5) in the middle
  double offset;    // offset for odd rows
  double rescale = 1.0 / (Nh * deltay);

  std::vector<double> points = {
      0, shift, 1, shift, 1, Nh * deltay + shift, 0, Nh * deltay + shift};
  std::vector<uint> segs = {0, 1, 1, 2, 2, 3, 3, 0};

  for (int i = 0; i <= Nh; ++i) // y coords
  {
    offset = (i % 2 != 0) ? deltax / 2 : 0; // shift at odd rows
    if (i % 2 != 0 && i < Nh) {             // additional point at odd rows
      points.push_back(0);
      points.push_back(i * deltay + shift);
    }
    for (int j = 0; j <= N; ++j) // x coords
    {
      if ((i == 0 && j == 0) || (i == Nh && j == N) ||
          (i == Nh && j + offset == 0) || (i == 0 && j == N))
        continue; // skip corners
      points.push_back(deltax * j +
                       (j == N ? 0 : offset)); // rightmost point within bound
      points.push_back(i * deltay + shift);
    }
  }

  for (int i = 1; i < points.size(); i += 2)
    points[i] *= rescale; // bring y in range [0,1]
  std::string flags;
  triangle_wrap(points, segs, {}, 0, flags.c_str(), coords_out, tris_out);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void anisotropic_tri(std::vector<double> &coords_out,
                     std::vector<uint> &tris_out, const uint N,
                     const double min_angle, const double y_anisotropy) {
  const float max_area = 2.1 * y_anisotropy / pow(N, 2);
  std::cout << y_anisotropy << " " << max_area << std::endl;
  std::string flags =
      // "Qq" + std::to_string(min_angle);
      "Qq" + std::to_string(min_angle) + "a" + std::to_string(max_area);


  // sampling
  const double& a = y_anisotropy;
  std::vector<double> points;
  std::vector<uint> edges;

  std::vector<double> xsamples = {0}, ysamples = {0};
  // const double l = sqrt(max_area * (4 / sqrt(3)));

  const double tolerance = .45;
  const double avgStep = a/static_cast<double>(N);

  // x in [0,a]
  double x = 0;
  while (x < a) {
    // random number in [-1,1]
    const double r = (static_cast<double>(rand()) /
      static_cast<double>(RAND_MAX)) * 2 - 1;
    x += avgStep * (1 + r * tolerance);
    if (x > a - avgStep*tolerance) x = a;
    xsamples.push_back(x);
  }

  // y in [0,1]
  double y = 0;
  while (y < 1) {
    // random number in [-1,1]
    const double r = (static_cast<double>(rand()) /
      static_cast<double>(RAND_MAX)) * 2 - 1;
    y += avgStep * (1 + r * tolerance);
    if (y > 1 - avgStep*tolerance) y = 1;
    ysamples.push_back(y);
  }

  const uint xN = xsamples.size(), yN = ysamples.size();
  // Generate points 
  for (int i = 0; i < xN; ++i) {
    points.push_back(0);
    points.push_back(xsamples.at(i));
  }
  for (int i = 0; i < yN; ++i) {
    points.push_back(ysamples.at(i));
    points.push_back(a);
  }
  for (int i = xN-1; i >= 0; --i) {
    points.push_back(1);
    points.push_back(xsamples.at(i));
  }
  for (int i = yN-1; i >= 0; --i) {
    points.push_back(ysamples.at(i));
    points.push_back(0);
  }
  // Generate edges
  edges.push_back(0);
  for (uint j = 1; j < 2*(xN+yN); ++j) {
    edges.push_back(j);
    edges.push_back(j);
  }
  edges.push_back(0);

  triangle_wrap(points, edges, {}, 0, flags.c_str(), coords_out,
                tris_out);
  for (int i = 1; i < coords_out.size(); i+=3)
    coords_out[i] /= y_anisotropy;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

double my_rand() {
  double r;
  do
    r = (double)rand() / RAND_MAX;
  while (r == 0.0 || r == 1.0);
  return r;
}

void non_unif_Poisson2D(int n, std::vector<double> &buf)
// non-uniform Poisson disc sampling of n points in [0,1]x[0,1]
// with non-uniform density inv. prop. to distance of point from origin
{
  double scale = n / 3;
  int edgesamples = 10 * (int)sqrt(n);
  int nsamples = 40000 - edgesamples;
  int count = 0;
  double x, y, r, rj, d;

  std::vector<double> bigbuf;
  std::vector<bool> flagbuf(nsamples, true);

  for (int i = 0; i < edgesamples; i++) {
    bigbuf.push_back(my_rand());
    bigbuf.push_back(0);
  }
  for (int i = 0; i < edgesamples; i++) {
    bigbuf.push_back(my_rand());
    bigbuf.push_back(1);
  }
  for (int i = 0; i < edgesamples; i++) {
    bigbuf.push_back(0);
    bigbuf.push_back(my_rand());
  }
  for (int i = 0; i < edgesamples; i++) {
    bigbuf.push_back(1);
    bigbuf.push_back(my_rand());
  }
  for (int i = 0; i < 2 * nsamples; i++)
    bigbuf.push_back(my_rand());

  for (int i = 0; i < bigbuf.size(); i += 2) {
    if (!flagbuf[i / 2])
      continue;
    x = bigbuf[i];
    y = bigbuf[i + 1];
    r = (x * x + y * y) / scale;
    for (int j = 0; j < bigbuf.size(); j += 2) {
      if (i == j)
        continue;
      if (!flagbuf[j / 2])
        continue;
      rj = (bigbuf[j] * bigbuf[j] + bigbuf[j + 1] * bigbuf[j + 1]) / scale;
      d = ((x - bigbuf[j]) * (x - bigbuf[j]) +
           (y - bigbuf[j + 1]) * (y - bigbuf[j + 1]));
      if (d < std::min(r, rj))
        flagbuf[j / 2] = false;
    }
  }

  count = 0;
  for (int i = 0; i < flagbuf.size(); i++) {
    if (flagbuf[i]) {
      buf.push_back(bigbuf[2 * i]);
      buf.push_back(bigbuf[2 * i + 1]);
      count++;
    }
    if (count >= n)
      break;
  }
  // std::cout << "Sampled " << count << " points out of desired " << n <<
  // std::endl;
}

void non_uniform_tri(std::vector<double> &coords_out,
                     std::vector<uint> &tris_out, int N) {
  double min_angle = 2.0; // this could be a parameter
  std::vector<double> points = {0, 0, 1, 0, 1, 1, 0, 1};
  std::vector<uint> segs = {0, 1, 1, 2, 2, 3, 3, 0};

  non_unif_Poisson2D(N * N, points);

  std::string flags; //= "Qq" + std::to_string(min_angle);
  triangle_wrap(points, segs, {}, 0, flags.c_str(), coords_out, tris_out);
}

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

void make_domain(std::vector<double> &coords_out, std::vector<uint> &tris_out,
                 const double max_area, const double min_angle = 20) {
  std::string flags =
      "Qq" + std::to_string(min_angle) + "a" + std::to_string(max_area);
  triangle_wrap({0, 0, 1, 0, 1, 1, 0, 1}, {0, 1, 1, 2, 2, 3, 3, 0}, {}, 0,
                flags.c_str(), coords_out, tris_out);
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

void twist_disk(double delta, std::vector<double> &coords_out)
{
  uint n=coords_out.size()/3;
  for (uint i=0;i<n;i++) {
    double x = coords_out[3*i], y = coords_out[3*i+1];
    double rho = sqrt(x*x+y*y);
    double theta = atan2(y,x)+delta*rho;
    coords_out[3*i] = rho * cos(theta);
    coords_out[3*i+1] = rho * sin(theta);
  }
}

void dilate_sin(double amp, std::vector<double> &coords_out)
{
  uint n=coords_out.size()/3;
  for (uint i=0;i<n;i++) {
    double x = coords_out[3*i], y = coords_out[3*i+1];
    double theta = atan2(y,x);
    double rho = sqrt(x*x+y*y);
    if (rho>0) rho += amp*(1-fabs(rho-0.5))*sin(18*theta);
    coords_out[3*i] = rho * cos(theta);
    coords_out[3*i+1] = rho * sin(theta);
  }
}



void write_soup(const char *filename,const vector<double> &coords, const vector<uint> &tris)
{
  ofstream f(filename);
  uint ntris = tris.size()/3;
  f << ntris << endl;
  for (int i=0;i<ntris;i++ ) {
    f << coords[3*tris[3*i]] << " " << coords[3*tris[3*i]+1] << " "
      << coords[3*tris[3*i+1]] << " " << coords[3*tris[3*i+1]+1] << " "
      << coords[3*tris[3*i+2]] << " " << coords[3*tris[3*i+2]+1] 
      << endl;
  }
  f.close();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void make_grid(std::vector<double> &coords_out, std::vector<uint> &tris_out,
               int N) {
  double delta = 1.0 / N;
  std::vector<double> points = {0, 0, 0, 1, 1, 1, 1, 0};
  std::vector<uint> segs = {0, 1, 1, 2, 2, 3, 3, 0};
  for (int i = 0; i <= N; ++i) {
    for (int j = 0; j <= N; ++j) {
      if ((i == 0 && j == 0) || (i == 0 && j == N) || (i == N && j == 0) ||
          (i == N && j == N)) {
        continue;
      }
      points.push_back(i * delta);
      points.push_back(j * delta);
    }
  }
  std::string flags;
  triangle_wrap(points, segs, {}, 0, flags.c_str(), coords_out, tris_out);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void make_triangulation(std::vector<double> &coords_out,
                        std::vector<uint> &tris_out, int mode, int N,
                        const double y_anisotropy) {

  const double min_angle = 0; // for Delaunay and anisotropic

  switch (mode) {
  case 0:
    regular_tri(coords_out, tris_out, N);
    break;
  case 1:
    make_domain(coords_out, tris_out, 1 / (1.6 * pow(N, 2)));
    break;
  case 2:
    non_uniform_tri(coords_out, tris_out, N);
    break;
  case 3:
    anisotropic_tri(coords_out, tris_out, N, min_angle, y_anisotropy);
    break;
  case 4:
    make_disk(N,coords_out,tris_out);
    break;
  }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                                 GUI utility

int main(int argc, char **argv) {
  std::vector<double> coords_out;
  std::vector<uint> tris_out;
  std::vector<uint> quads{};
  uint N = 10;  // default
  double delta = 10;  // rotations per step
  uint n_disks = 1;   // # steps-files
  double amp = 0;     // amplification for radial displacement
  if (argc > 1) N = std::stoi(argv[1]);
  if (argc > 2) n_disks = std::stoi(argv[2]);
  if (argc > 3) delta = std::stod(argv[3]);
  if (argc > 4) amp = std::stod(argv[4]);
  make_triangulation(coords_out, tris_out, 4, N, delta);
  std::string namebase = "../output/disk" + std::to_string(N);
  std::string name = namebase  + ".obj";
  write_OBJ(name.c_str(), coords_out, tris_out, quads);
  name = namebase  + ".msh";
  write_soup(name.c_str(), coords_out, tris_out);
  for (uint i=0;i<n_disks;i++) {
    twist_disk(delta*M_PI/180.0,coords_out);
    if (amp>0) dilate_sin(amp,coords_out);
    name = namebase + "_t" + std::to_string(static_cast<int>(delta*(i+1))) + ".obj";
    write_OBJ(name.c_str(), coords_out, tris_out, quads);
    name = namebase + "_t" + std::to_string(static_cast<int>(delta*(i+1))) + ".msh";
    write_soup(name.c_str(), coords_out, tris_out);
  }
}
