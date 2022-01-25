#ifndef DUNE_ULTRAWEAK_FULL_BLOCK_INVERT_HH
#define DUNE_ULTRAWEAK_FULL_BLOCK_INVERT_HH

template<typename M, typename X>
class FullBlockInvert : public Dune::InverseOperator<X,X> {
public:

  FullBlockInvert(M& mat_) : mat(mat_) {
    for(auto& block : mat)
      block.begin()->invert();
  }

  void apply(X& x, X& b, Dune::InverseOperatorResult& res) {
    mat.mv(b, x);
    res.iterations = 1;
    res.converged = true;
  }

  void apply(X& x, X& b, double reduction, Dune::InverseOperatorResult& res) {
    this->apply(x, b, res);
  }

  Dune::SolverCategory::Category category() const {
    return Dune::SolverCategory::sequential;
  }

  const M& getmat() const {
    return mat;
  }

private:
  M& mat;
};

#endif // DUNE_ULTRAWEAK_FULL_BLOCK_INVERT_HH
