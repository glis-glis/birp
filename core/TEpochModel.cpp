
//
// Created by Liam on 23/08/22.
//

#include "TEpochModel.h"
#include "TBirpPrior.h"

using namespace coretools::instances;

TEpochModel::TEpochModel(stattools::TParameterTyped<TypeGamma, 1> *Gamma, stattools::TParameterTyped<TypeTime, 1>* TimesOfChange,
                        const std::string &Filename)
        : _gamma(Gamma), _timesOfChange(TimesOfChange), _filename(Filename) {



}

void TEpochModel::initialize() {
    if (!_timesOfChange->hasFixedInitialValue()){  // if user did not set initial values
		_numEpochs = 1;
	}
	_numEpochs = _timesOfChange->totalSize() + 1;
	_gamma->initStorage({_numEpochs}, {std::make_shared<coretools::TNamesIndices>()});
}


std::vector<size_t> TEpochModel::_getEpochStarts(const TUniqueContainer<size_t> &Timepoints) const {
    assert(!Timepoints.empty());
    std::vector<size_t> epochStarts = {0};
    if (_timesOfChange->totalSize()!=0) {
        size_t maxTimeOfChange = _timesOfChange->max_value();
        size_t maxTime         = std::max(Timepoints.back(), maxTimeOfChange + 1);
		for (size_t i = 0; i < _timesOfChange->totalSize(); i++){
			epochStarts.push_back(_timesOfChange->value(i));
		}
        epochStarts.push_back(maxTime);
    } else {
        epochStarts.push_back(*max_element(Timepoints.cbegin(), Timepoints.cend()));
    }
    return epochStarts;
}

// static void _readTimesOfChange(){}


void TEpochModel::_fillRho(const TUniqueContainer<size_t> &Timepoints) {
    std::vector<size_t> epochStarts = _getEpochStarts(Timepoints);
    _rho.zeros((int)Timepoints.size(), (int)_numEpochs);
    // fill rhos
    for (size_t t = 0; t < Timepoints.size(); ++t) {
        for (size_t e = 1; e < epochStarts.size(); ++e) {
            if (Timepoints[t] <= epochStarts[e - 1]) {
                _rho(t, e - 1) = 0;
            } else if (epochStarts[e - 1] < Timepoints[t] && Timepoints[t] < epochStarts[e]) {
                _rho(t, e - 1) = Timepoints[t] - epochStarts[e - 1];
            } else {
                _rho(t, e - 1) = epochStarts[e] - epochStarts[e - 1];
            }
        }
    }
}


double TEpochModel::_calculateLogPhi_k(size_t k) const {
    double sum = 0.0;
    for (size_t m = 0; m < _numEpochs; m++) { sum += (double)_rho(k, m) * _gamma->value(m); }
    return sum;
}

std::vector<double> TEpochModel::_calculatePhi(size_t numTimepoints) const {
    std::vector<double> phi(numTimepoints);
    for (size_t k = 0; k < numTimepoints; k++) { phi[k] = _calculateLogPhi_k(k); }
    coretools::standardizeZeroMean(phi);
    for (auto &p : phi) { p = exp(p); } // make logPhi to phi
    return phi;
}

/*
coretools::TMatrix<double> TEpochModel::_calculateP(size_t numTimepoints, size_t numMethods, const TSpecies &Data) const {
   coretools::TMatrix<double> p((int)numMethods, (int)numTimepoints);
   for (size_t i = 0; i < numMethods; i++) {
       double sum = 0.0;
       for (size_t k = 0; k < numTimepoints; k++) {
           // calculate sigma and fill into p (will normalize later)
           const double sigma = _phi[k] * Data() (i,k)   _calculateSigma_ik(i, k, _phi[k]);
           p(i, k)            = sigma;
           sum += sigma;
       }
       // now normalize by row sum!
       for (size_t k = 0; k < numTimepoints; k++) { p(i, k) /= sum; }
   }
   return p;
}

void TEpochModel::_evaluateUpdateGamma(size_t m) {
  // calculate probabilities p for new gamma
  coretools::TMatrix<double> p; // = _calculateP();


double calculateLL(const TDataSeries &Data) {

}


const double newLL  = calculateLL(p);
const double LLRatio = newLL - _oldLL;

const double newLogJeffPrior = _calculateLogJeffreysPrior(p);
const double logPriorRatio   = newLogJeffPrior - _oldLogJeffPrior;

const double logH = LLRatio + logPriorRatio;
if (_gamma->acceptOrReject(m, logH)) {
  // accepted -> update LL and jeffrey's Prior
  _oldLL           = newLL;
  _oldLogJeffPrior = newLogJeffPrior;
}
}

void TEpochModel::updateGamma(){
for (size_t m = 0; m < _numEpochs; m++) {
   _gamma->updateSpecificIndex(m);
   _evaluateUpdateGamma(m);
}
}

*/

// getters for test functions
size_t TEpochModel::rho(size_t row, size_t col) const {
   return _rho(row,col);
}






/*
namespace
/*anonymous*/ /*
  {
     class Sums{
     private:
         double LL{};
         TypeCounts sumN{};
         double sumSigma{};
         TypeIndex location{};
         TypeIndex method{};

         void _setIndex(const TDataPoint & D){
             location = D.location;
             method = D.method;
         }

         void _initSums(const TDataPoint& D, const std::vector<double> & Phi){
             double sigma = Phi[D.timepoint] * D.effort;
             sumSigma = sigma;
             sumN = D.counts;
             LL += D.counts * log(sigma);
         }

         void _addToSums(const TDataPoint& D, const std::vector<double> & Phi){
             double sigma = Phi[D.timepoint] * D.effort;
             sumSigma += sigma;
             sumN += D.counts;
             LL += D.counts * log(sigma);
             _setIndex(D);
         }

         void _normalize(const TDataPoint& D, const std::vector<double> & Phi){
             LL -= sumN * log(sumSigma);
             _initSums(D, Phi);
         }

     public:
         Sums(const TDataPoint& D, const std::vector<double> & Phi){
             LL = 0.0;
             _initSums(D, Phi);
             _setIndex(D);
         }

         void add(const TDataPoint& D, const std::vector<double> & Phi){
             if(D.location != location || D.method != method){
                 _normalize(D, Phi);
             } else {
                 _addToSums(D, Phi);
             }
         }

         double getLL(){
             return LL;
         }
     };
 }

 double TEpochModel::calculateLL(const TData & Data){
     auto d = Data.cbegin();
     Sums sums(*d, _phi);
     for(; d != Data.cend(); ++d){
         sums.add(*d, _phi);
     }
     return sums.getLL();
 }


 */
