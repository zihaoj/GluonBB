#ifndef Helpers_Helpers_H
#define Helpers_Helpers_H

#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODTracking/Vertex.h"
#include <functional>
#include <cmath>
#include "xAODBTagging/BTagging.h"
#include "xAODJet/Jet.h"
#include "xAODJet/JetAttributes.h"


float MV2(const xAOD::Jet* jet, std::string flavor);
float BTAGIP(const xAOD::Jet* inputjet, std::string flavor);
int   MV2_benchmark(const xAOD::Jet* jet, std::string flavor, int name);
float dphi(float phi1, float phi2);
float dR(float eta1, float phi1, float eta2, float phi2);
float dR_vertices(float x1, float y1, float z1, float x2, float y2, float z2);
bool came_from_b_quark(const xAOD::TruthParticle* p);
bool has_b_hadron_child(const xAOD::TruthParticle* p);
template <class T> struct sort_pt : std::binary_function <T,T,bool> {
  bool operator() (T* x, T* y) const {return x->pt()>y->pt();}
};
bool svtx_inside_jet(double x, double y, double z, const xAOD::Vertex* pv, const xAOD::Jet* jet, double dR_thresh);
//int count_vertices(std::vector<double> x_vert, std::vector<double> y_vert, std::vector<double> z_vert, const xAOD::Vertex* pv, const xAOD::Jet* jet, double dR_thresh = 0.2);
int count_vertices(std::vector<double> x_vert, std::vector<double> y_vert, std::vector<double> z_vert, const xAOD::Jet* jet, double dR_thresh = 0.2);
int count_bquarks(std::vector<double> eta, std::vector<double> phi, std::vector<int> pdgid, const xAOD::Jet* jet, double dR_thresh = 0.2);


//inline double d0UncertaintyBeamSpot2(double track_phi0, double beam_sigma_x, double beam_sigma_y, double beam_sigma_xy){
//  double sin_phi = sin(track_phi0);
//  double cos_phi = cos(track_phi0);
//  double d0_uncert2 = sin_phi *(  sin_phi*pow(beam_sigma_x,2)
//				  -cos_phi*beam_sigma_xy)
//    +cos_phi *( cos_phi *pow(beam_sigma_y,2) - sin_phi*beam_sigma_xy);
//  return d0_uncert2;
//    
//}
//
//double d0significance(const xAOD::TrackParticle *tp, double beam_sigma_x, double beam_sigma_y, double beam_sigma_xy) {
//   //checkTPAndDefiningParamCov(tp);
//   double d0 = tp->d0();
//   // elements in definingParametersCovMatrixVec should be : sigma_d0^2, sigma_d0_z0, sigma_z0^2
//   double sigma_d0 = tp->definingParametersCovMatrixVec().at(0);
//   return d0/sqrt(sigma_d0+d0UncertaintyBeamSpot2(tp->phi(),beam_sigma_x, beam_sigma_y, beam_sigma_xy));
//}

#endif
