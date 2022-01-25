// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_ULTRAWEAK_CONSTRAINTS_HH
#define DUNE_ULTRAWEAK_CONSTRAINTS_HH

template<typename Mat, typename Rhs, typename CC>
void eliminate_cols(Mat& mat, Rhs& rhs, const CC& cc) {
  using BlockType = typename Mat::row_type::member_type;

  for (const auto& dof : cc) {
    auto idx = dof.first[0]; // TODO: not working for true multiindices!
    for (auto& row : mat) { // TODO: separate those loops
      auto it = row.find(idx);
      if(!it.equals(row.end())) {
        // TODO: save data
        *it = BlockType(0.0);
        rhs[it.index()] = 0.0;
      }
    }
  }
}

// this is some hardcoded stuff that will likely not work in other contexts
template<typename Mat, typename CC>
void eliminateAmatCols(Mat& aMat, const CC& cc) {
  using BlockType = typename Mat::row_type::member_type;

  for (const auto& dof : cc) {
    auto idx = dof.first[0]; // TODO: not working for true multiindices!
    for (auto row=aMat.begin(); row!=aMat.end(); row++) {
      if(row.index() != idx) {
        auto it = row->find(idx);
        if(!it.equals(row->end()))
          *it = BlockType(0.0);
      }
    }
  }
}

#endif  // DUNE_ULTRAWEAK_CONSTRAINTS_HH
