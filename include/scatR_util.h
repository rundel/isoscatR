#ifndef _SCATR_UTIL_H
#define _SCATR_UTIL_H

#include <RcppArmadillo.h>
#include "scatR_structs.h"

void open_allelefiles(GlobalParams &p, GlobalOptions &opt);
void close_allelefiles(GlobalParams &p, GlobalOptions &opt);

void open_cvfiles(GlobalParams &p, GlobalOptions &opt);
void close_cvfiles(GlobalParams &p, GlobalOptions &opt);

double calc_multinomial_loglik( const arma::mat &theta, const arma::mat &count, const arma::colvec &sumCount);
double calc_multivar_normal_loglik(const arma::mat &x, const arma::mat &mu, const arma::mat &sigmaInv, double sigmaDet);

arma::mat calc_f(const arma::mat &theta);

arma::mat calc_rotation_mat(double angle);
arma::mat calc_stretch_mat(double ratio);

void outputAccepts(GlobalParams &p, GlobalOptions &opt);
void outputTuning(GlobalParams &p, GlobalOptions &opt);
void outputPerformance(GlobalParams &p, GlobalOptions &opt);

double accept_ratio(unsigned int accept, unsigned int attempt);

double great_circle_dist(double x1, double y1, double x2, double y2);
arma::mat distance_mat(arma::mat const &m);
arma::mat distance_mat(arma::colvec const &xc, arma::colvec const &yc);

double isLeft(double x0, double y0, double x1, double y1,  double x2, double y2);
int isInsideBoundary( double x, double y, arma::mat boundary);

#endif