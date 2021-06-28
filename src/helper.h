#ifndef HELPER_H
#define HELPER_H

template <typename T = double> bool areEqual(T a, T b) {
  return std::fabs(a - b) < 1e-6;
}

template <typename T = double>
void print_vec(const vector<T> &v, const string prefix = "") {
  cout << prefix << " SIZE=" << v.size() << " : ";
  for (auto x : v)
    cout << x << ',';
  cout << endl;
}

template<typename T = double>
void print_arr( const double* a, const size_t n, const string prefix="" ) {
    cout << prefix << " SIZE=" << n << " : ";
    for (size_t i = 0; i < n; i++) 
        cout << a[i] << ',';
    cout << endl;
}

#endif /* end of include guard: HELPER_H */
