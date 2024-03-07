/**
 * @file    sparse_vector_test.cc
 * @brief   Implementation tests for a sparse_vector object.
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    April 22, 2014
 * @version 1.0
 *
 * The test program is made of some infrastucture (mostly the TestManagerClass
 * class) and some test operations.
 *
 * The test operations are derived by actions::BaseAction class; while they are
 * not documented, their content is a good showcase of features of the
 * sparse_vector class.
 */

// C/C++ standard libraries
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept> // std::out_of_range
#include <string>
#include <utility> // std::make_pair()

// LArSoft (larcore) libraries
#include "lardataobj/Utilities/sparse_vector.h"

template <typename T>
T abs(T value)
{
  return (value < T(0)) ? -value : value;
}

/// Prints a container into a stream
template <typename CONT>
std::ostream& PrintVector(const CONT& sv, std::ostream& out = std::cout)
{
  out << "(" << sv.size() << ") {";
  for (auto value : sv)
    out << " " << value;
  return out << " }" << std::flush;
} // PrintVector()

/// Converts a container into a string
template <typename CONT>
std::string VectorToString(const CONT& sv)
{
  std::ostringstream sstr;
  PrintVector(sv, sstr);
  return sstr.str();
} // VectorToString()

/// Prints a sparse vector highlighting the ranges
template <typename T>
std::ostream& PrintVectorRanges(const lar::sparse_vector<T>& sv, std::ostream& out = std::cout)
{
  out << "(" << sv.size() << ") [" << sv.get_ranges().size() << "] {";
  size_t i = 0;
  for (auto range : sv.get_ranges()) {
    while (i++ < range.begin_index())
      out << " " << lar::sparse_vector<T>::value_zero;
    out << " [";
    for (auto value : range)
      out << " " << value;
    out << " ]";
    i += range.size() - 1;
  }
  while (i++ < sv.size())
    out << " " << lar::sparse_vector<T>::value_zero;
  return out << " }" << std::flush;
} // PrintVectorRanges()

//------------------------------------------------------------------------------
namespace actions {
  template <typename T>
  class BaseAction;
} // namespace actions

/// Test manager class
template <typename T>
class TestManagerClass {
public:
  typedef T Data_t;
  typedef TestManagerClass<Data_t> This_t;
  typedef std::vector<Data_t> Vector_t;
  typedef lar::sparse_vector<Data_t> SparseVector_t;
  typedef actions::BaseAction<Data_t> Action_t;

  /// Data structure to document failure of a single action
  struct FailureInfo_t {
    unsigned int nAction{0}; ///< ID of the action
    std::string description; ///< description of the failed action
    unsigned int nErrors{0}; ///< number of detected errors

    FailureInfo_t() {}
    FailureInfo_t(This_t& tc, const Action_t& action, unsigned int errors)
      : nAction(tc.nAction), description(action.description(tc)), nErrors(errors)
    {}
  }; // FailureInfo_t

  Vector_t v{};                        ///< vector to operate on
  SparseVector_t sv{};                 ///< sparse vector to operate on
  unsigned int nAction{0};             ///< number of actions taken so far
  unsigned int nErrors{0};             ///< number of errors accumulated
  std::vector<FailureInfo_t> failures; ///< documentation of the failures

  int quietness{0}; ///< verbosity level (larger is less verbose)

  std::ostream& out; ///< output stream

  /// Constructor, specify an output stream
  TestManagerClass(std::ostream& output_stream = std::cout) : out(output_stream) {}

  /// Returns the current size of the STL vector
  size_t current_vector_size() const { return v.size(); }

  /// Checks if the sparse vector is valid
  bool validate() const
  {
    const auto ranges = sv.get_ranges();
    // a sparse vector with no non-null elements can't be detected invalid
    if (ranges.empty()) {
      //  std::cout << "[vector is void]" << std::endl;
      return true;
    }
    int i = 0;
    auto iNext = ranges.cbegin(), rend = ranges.cend();
    while (iNext != rend) {
      auto iRange = iNext++;
      if (iRange->empty()) {
        out << "[range " << i << " is empty]" << std::endl;
        return false;
      }
      if (iNext != rend) {
        if (!(*iRange < *iNext)) {
          out << "[range " << i << " comes after range " << (i + 1) << "]" << std::endl;
          return false;
        }
        if (!iRange->separate(*iNext)) {
          out << "[range " << i << " is overlapping or adjacent to " << (i + 1) << "]" << std::endl;
          return false;
        }
      }
      ++i;
    } // while
    if (sv.size() < ranges.back().end_index()) {
      out << "[size is " << sv.size() << ", it should be at least " << ranges.back().end_index()
          << "]" << std::endl;
      return false;
    }
    return true;
  } // validate()

  bool match() const { return std::equal(v.begin(), v.end(), sv.begin()); }

  unsigned int check()
  {
    unsigned int local_errors = 0;
    if (!validate()) ++local_errors;
    if (!match()) {
      out << "*** Mismatch: " << std::string(58, '*') << std::endl;
      out << "vector:        " << VectorToString(v) << std::endl;
      out << "sparse vector: " << VectorToString(sv) << std::endl;
      out << std::string(72, '*') << std::endl;
      ++local_errors;
    }
    return local_errors;
  } // check()

  void run(const Action_t& action)
  {
    action(v);
    action(sv);
  } // run()

  /**
   * @brief Perform the specified action and evaluate the result
   * @param action the action to be performed
   * @param autoCorrect if on failure the sparse vector is copied from vector
   * @return number of detected errors
   */
  unsigned int perform(const Action_t& action, bool autoCorrect = true)
  {
    out << "[" << ++nAction << "] ";
    action.describe(*this);
    out << std::endl;
    run(action);
    if (quietness <= 0) {
      out << "    => ";
      PrintVectorRanges(sv) << std::endl;
    }
    unsigned int errors = check();
    if (action.expected_errors != 0) {
      out << "  (" << action.expected_errors << " error(s) expected)" << std::endl;
    }
    if ((int)errors != action.expected_errors) {
      int actual_errors = abs(errors - action.expected_errors);
      nErrors += actual_errors;
      failures.emplace_back(*this, action, actual_errors);
      if (autoCorrect) {
        recover();
        if (!match()) out << "   *** recovery from error failed ***" << std::endl;
      } // if autoCorrect
    }   // if errors
    return errors;
  } // perform()

  int quiet(int nq = 0)
  {
    int q = quietness;
    quietness = nq;
    return q;
  }

  unsigned int operator()(const Action_t& action) { return perform(action); }

  /// Set the sparse vector as a copy of the vector (lossy!)
  void recover()
  {
    sv.assign(v);
    sv.optimize();
  }

  /// Set the vector as a copy of the sparse vector
  void mirror() { v.assign(sv.cbegin(), sv.cend()); }

  /// Print a summary of the failures
  int summary() const
  {
    out << nAction << " tests performed";
    if (nErrors == 0) { out << " successfully." << std::endl; }
    else {
      out << ", " << nErrors << " unexpected errors occurred in " << failures.size()
          << " of them:" << std::endl;
      for (const FailureInfo_t& info : failures) {
        out << "  [" << info.nAction << "] " << info.description << " => " << info.nErrors
            << " errors" << std::endl;
      } // for
    }   // if errors
    return nErrors;
  } // summary()

}; // TestManagerClass<>

//------------------------------------------------------------------------------
/// Namespace containing all the test actions.
namespace actions {

  /**
   * @brief Base class for testing actions
   * @param T the type stored in the containers
   *
   * An action class performs a definite action on both a STL vector (Vector_t)
   * and on a sparse_vector (SparseVector_t).
   * The interface includes:
   * - a describe() function describiing the detail of the action
   * - a operator(SparseVector_t&) on a sparse vector
   * - a operator(Vector_t&) on a STL vector
   * In most of the cases, the action on the STL vector mirrors the one in
   * the sparse one. Therefore each action shows how to do the same thing
   * in the two different containers.
   */
  template <typename T>
  class BaseAction {
  public:
    using This_t = BaseAction<T>;
    typedef TestManagerClass<T> TestClass_t;
    typedef typename TestClass_t::Data_t Data_t;     ///< type stored
    typedef typename TestClass_t::Vector_t Vector_t; ///< type of STL vector
    typedef typename TestClass_t::SparseVector_t SparseVector_t;
    ///< type of sparse vector

    int expected_errors{0}; ///< number of errors expected for this action

    BaseAction() {}
    BaseAction(int exp_err) : expected_errors(exp_err) {}

    virtual ~BaseAction() {}

    This_t& set_expected_errors(int exp_err)
    {
      expected_errors = exp_err;
      return *this;
    }

    /// Print a description of this action into the specified steam.
    void describe(TestClass_t& tc, std::ostream& out) const { doDescribe(tc, out); }

    /// Print a description of this action in the test class stream.
    void describe(TestClass_t& tc) const { return doDescribe(tc, tc.out); }

    /// Returns a string with a description of this action.
    std::string description(TestClass_t& tc) const
    {
      std::ostringstream sstr;
      describe(tc, sstr);
      return sstr.str();
    }

    /// Action performed on a STL vector.
    void operator()(Vector_t& v) const { actionOnVector(v); }

    /// Action performed on a sparse vector.
    void operator()(SparseVector_t& v) const { actionOnSparseVector(v); }

  protected:
    using v_citer_t = typename Vector_t::const_iterator;

    virtual void doDescribe(TestClass_t&, std::ostream& out) const { out << "no action"; }

    virtual void actionOnVector(Vector_t&) const {}
    virtual void actionOnSparseVector(SparseVector_t&) const {}

    template <typename ITER>
    static ITER skipRange(ITER start, ITER const end);
    template <typename ITER>
    static ITER skipVoid(ITER start, ITER const end);
    template <typename ITER>
    static ITER rangeStart(ITER const begin, ITER last);
    template <typename ITER>
    static ITER voidStart(ITER const begin, ITER last);
    template <typename Vector>
    static auto findRangeBorders(Vector& v, std::size_t pos)
    {
      if ((pos >= v.size()) || isVoid(v[pos])) return std::make_pair(v.end(), v.end());
      auto rbegin = rangeStart(v.begin(), v.begin() + pos);
      auto rend = skipRange(rbegin, v.end());
      return std::make_pair(rbegin, rend);
    }

    /// Returns iterator to the first element of the i-th range after start.
    template <typename ITER>
    static ITER findRangeStart(Vector_t const& v, std::size_t i, ITER start);
    template <typename ITER>
    static ITER findRangeStart(Vector_t const& v, std::size_t i)
    {
      return findRangeStart(v, i, v.cbegin());
    }

    /// Returns iterator to the first void element of the i-th void area
    /// after start.
    template <typename ITER>
    static ITER findVoidStart(Vector_t const& v, std::size_t i, ITER start);

    static bool isVoid(Data_t value) { return SparseVector_t::is_zero(value); }

  }; // BaseAction<>

  template <typename T>
  class Print : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

  protected:
    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "print the vectors";
    }

    virtual void actionOnVector(Vector_t& v) const override { doPrintVector(v); }
    virtual void actionOnSparseVector(SparseVector_t& v) const override { doPrintVector(v); }

    template <typename Vector>
    static void doPrintVector(Vector const& v)
    {
      PrintVector(v) << std::endl;
    }

  }; // Print<>

  template <typename T>
  class PrintSparseVector : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

  protected:
    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "print the sparse vector";
    }

    virtual void actionOnSparseVector(SparseVector_t& v) const override
    {
      PrintVector(v) << std::endl;
    }

  }; // PrintSparseVector<>

  template <typename T>
  class PrintRange : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

    size_t position;

    PrintRange(size_t pos) : position(pos) {}

  protected:
    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "print the range at position " << position;
    }

    virtual void actionOnSparseVector(SparseVector_t& v) const override
    {
      try {
        std::string range_str = VectorToString(v.find_range(position));
        std::cout << "   position " << position << " (value: " << v[position] << " is within range "
                  << range_str << std::endl;
      }
      catch (std::out_of_range const&) {
        std::cout << "   no range at position " << position << std::endl;
      }
    }

  }; // PrintRange<>

  template <typename T>
  class PrintNonVoid : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

  protected:
    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "print the elements which are not in the void";
    }

    virtual void actionOnVector(Vector_t& v) const override
    {
      typedef typename Base_t::SparseVector_t::value_type (*cmp_t)(
        typename Base_t::SparseVector_t::value_type);
      std::cout << "Non-zero elements in vector:        "
                << (v.size() -
                    std::count_if(v.begin(), v.end(), cmp_t(&Base_t::SparseVector_t::is_zero)))
                << std::endl;
    }
    virtual void actionOnSparseVector(SparseVector_t& v) const override
    {
      std::cout << "Non-void elements in sparse vector: " << v.count() << std::endl;
    }

  }; // PrintNonVoid<>

  template <typename T>
  class Clear : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

  protected:
    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "clear the vectors";
    }

    virtual void actionOnVector(Vector_t& v) const override { doClearVector(v); }
    virtual void actionOnSparseVector(SparseVector_t& v) const override { doClearVector(v); }

    template <typename Vector>
    static void doClearVector(Vector& v)
    {
      v.clear();
    }

  }; // Clear<>

  template <typename T>
  class Resize : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

    size_t size;

    Resize(size_t new_size = 0) : size(new_size) {}

  protected:
    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "resize to " << size << " elements";
    }

    virtual void actionOnVector(Vector_t& v) const override { doResizeVector(v); }
    virtual void actionOnSparseVector(SparseVector_t& v) const override { doResizeVector(v); }

    template <typename Vector>
    void doResizeVector(Vector& v) const
    {
      v.resize(size);
    }

  }; // Resize<>

  template <typename T>
  class Truncate : public Resize<T> {
  public:
    using Base_t = Resize<T>;
    using typename Base_t::TestClass_t;

    Truncate(size_t new_size = 0) : Base_t(new_size) {}

  protected:
    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "truncate to " << Base_t::size << " elements";
    }

  }; // Truncate<>

  template <typename T>
  class ResizeWith : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

    size_t size;
    Data_t value;

    ResizeWith(size_t new_size, Data_t val) : size(new_size), value(val) {}

  protected:
    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "resize to " << size << " elements, filling with " << value;
    }

    virtual void actionOnVector(Vector_t& v) const override { doResizeVector(v); }
    virtual void actionOnSparseVector(SparseVector_t& v) const override { doResizeVector(v); }

    template <typename Vector>
    void doResizeVector(Vector& v) const
    {
      v.resize(size, value);
    }

  }; // ResizeWith<>

  template <typename T>
  class Assign : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

    Vector_t data;

    Assign(Vector_t new_data) : data(new_data) {}

  protected:
    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "assign ";
      PrintVector(data, out);
    } // describe()

    virtual void actionOnVector(Vector_t& v) const override { doAssignVector(v); }
    virtual void actionOnSparseVector(SparseVector_t& v) const override { doAssignVector(v); }

    template <typename Vector>
    void doAssignVector(Vector& v) const
    {
      v.assign(data.begin(), data.end());
    }

  }; // Assign<>

  template <typename T>
  class AssignMove : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

    Vector_t data;

    AssignMove(Vector_t new_data) : data(new_data) {}

  protected:
    virtual void actionOnVector(Vector_t& v) const override
    {
      Vector_t local(data);
      v.swap(local);
    }
    virtual void actionOnSparseVector(SparseVector_t& v) const override
    {
      SparseVector_t local(data);
      v.assign(std::move(local));
    }

    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "move ";
      PrintVector(data, out);
      out << " into the vector";
    } // describe()

  }; // AssignMove<>

  template <typename T>
  class Insert : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

    size_t position;
    Vector_t data;

    Insert(size_t pos, Vector_t new_data) : position(pos), data(new_data) {}

  protected:
    virtual void actionOnVector(Vector_t& v) const override
    {
      size_t max_size = std::max(v.size(), position + data.size());
      v.resize(max_size, 0);
      std::copy(data.begin(), data.end(), v.begin() + position);
    }
    virtual void actionOnSparseVector(SparseVector_t& v) const override
    {
      v.add_range(position, data);
    }

    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "add data vector ";
      PrintVector(data, out);
      out << " starting at position " << position;
    } // describe()

  }; // Insert<>

  template <typename T>
  class Add : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

    size_t position;
    Vector_t data;
    Data_t baseline = Data_t(0);

    Add(size_t pos, Vector_t new_data, Data_t baseline = Data_t(0))
      : position(pos), data(new_data), baseline(baseline)
    {}

  protected:
    virtual void actionOnVector(Vector_t& v) const override
    {
      size_t max_size = std::max(v.size(), position + data.size());
      v.resize(max_size, 0);
      // buggy: we have no way to recognise a void element from a 0 one
      std::transform(data.cbegin(),
                     data.cend(),
                     v.cbegin() + position,
                     v.begin() + position,
                     [this](auto a, auto b) { return a + ((b == 0) ? baseline : b); });
    }
    virtual void actionOnSparseVector(SparseVector_t& v) const override
    {
      v.combine_range(position, data, std::plus<Data_t>(), baseline);
    }

    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "increment by data in ";
      PrintVector(data, out);
      out << " with baseline " << baseline << " starting at position " << position;
    } // describe()

  }; // Add<>

  template <typename T>
  class Erase : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

    size_t first, last;

    Erase(size_t from, size_t to) : first(from), last(to) {}

  protected:
    virtual void actionOnVector(Vector_t& v) const override
    {
      std::fill(v.begin() + std::min(first, v.size()),
                v.begin() + std::min(last, v.size()),
                Base_t::SparseVector_t::value_zero);
    }
    virtual void actionOnSparseVector(SparseVector_t& v) const override
    {
      v.make_void(v.begin() + first, v.begin() + last);
    }

    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "render elements from " << first << " to " << last << " void";
    } // describe()

  }; // Erase<>

  template <typename T>
  class EraseRangeAt : public BaseAction<T> {
    /*
     * Note: `EraseRangeAt<>` prediction might be less than perfect since
     * std::vector does not have information about the ranges.
     * In that case, `set_expected_errors(1)` should be called before testing,
     * and `mirror()` after the test.
     */
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

    size_t position;

    EraseRangeAt(size_t pos) : position(pos) {}

  protected:
    virtual void actionOnVector(Vector_t& v) const override
    {
      auto&& [rstart, rend] = Base_t::findRangeBorders(v, position);
      std::fill(rstart, rend, SparseVector_t::value_zero);
    }

    virtual void actionOnSparseVector(SparseVector_t& v) const override
    {
      v.make_void_around(position);
    }

    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "void range containing position " << position;
    }

  }; // EraseRangeAt<>

  template <typename T>
  class Scale : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

    size_t position;
    Data_t factor;

    Scale(size_t pos, Data_t factor) : position(pos), factor(factor) {}

  protected:
    virtual void actionOnVector(Vector_t& v) const override
    {
      auto&& [rstart, rend] = Base_t::findRangeBorders(v, position);
      while (rstart < rend)
        *(rstart++) *= factor;
    }
    virtual void actionOnSparseVector(SparseVector_t& v) const override
    {
      auto i = v.find_range_number(position);
      if (i >= v.n_ranges()) return;
      for (auto& value : v.range_data(i))
        value *= factor;
    }

    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "scale data by a factor " << factor << " starting at position " << position;
    } // describe()

  }; // Scale<>

  // fun with C++17:
  template <typename T>
  explicit Scale(size_t pos, T factor) -> Scale<T>;

  template <typename T>
  class ScaleAll : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

    Data_t factor;

    ScaleAll(Data_t factor) : factor(factor) {}

  protected:
    virtual void actionOnVector(Vector_t& v) const override
    {
      for (auto& value : v)
        value *= factor;
    }
    virtual void actionOnSparseVector(SparseVector_t& v) const override
    {
      for (auto& r : v.iterate_ranges()) {
        for (auto& value : r)
          value *= factor;
      }
    }

    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "scale all data by a factor " << factor;
    }

  }; // ScaleAll<>

  template <typename T>
  explicit ScaleAll(T factor) -> ScaleAll<T>;

  template <typename T>
  class SetElement : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

    size_t position;
    Data_t value;

    SetElement(size_t pos, Data_t new_value) : position(pos), value(new_value) {}

  protected:
    virtual void actionOnVector(Vector_t& v) const override { v[position] = value; }
    virtual void actionOnSparseVector(SparseVector_t& v) const override
    {
      v.set_at(position, value);
    }

    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "set the element " << position << " to " << value;
    }

  }; // SetElement<>

  template <typename T>
  class UnsetElement : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

    size_t position;

    UnsetElement(size_t pos) : position(pos) {}

  protected:
    virtual void actionOnVector(Vector_t& v) const override
    {
      v[position] = Base_t::SparseVector_t::value_zero;
    }
    virtual void actionOnSparseVector(SparseVector_t& v) const override { v.unset_at(position); }

    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "turn the element " << position << " into the void";
    }

  }; // SetElement<>

  template <typename T>
  class PushBack : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

    Data_t value;

    PushBack(Data_t val) : value(val) {}

  protected:
    virtual void actionOnVector(Vector_t& v) const override { v.push_back(value); }
    virtual void actionOnSparseVector(SparseVector_t& v) const override
    {
      v.push_back(value, 0. /* SparseVector_t::value_zero */);
    }

    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "add element " << value << " (if non-zero)";
    }

  }; // PushBack<>

  template <typename T>
  class SetValue : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

    size_t position;
    Data_t value;

    SetValue(size_t pos, Data_t val) : position(pos), value(val) {}

  protected:
    virtual void actionOnVector(Vector_t& v) const override { doAssignValue(v); }
    virtual void actionOnSparseVector(SparseVector_t& v) const override { doAssignValue(v); }

    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "set element " << position << " to " << value;
    }

    template <typename Vector>
    void doAssignValue(Vector& v) const
    {
      v[position] = value;
    }

  }; // SetValue<>

  template <typename T>
  class FlipSign : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

    size_t position;

    FlipSign(size_t pos) : position(pos) {}

  protected:
    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "change the sign of element " << position;
    }

    virtual void actionOnVector(Vector_t& v) const override { doFlipValue(v); }
    virtual void actionOnSparseVector(SparseVector_t& v) const override { doFlipValue(v); }

    template <typename Vector>
    void doFlipValue(Vector& v) const
    {
      v[position] = -v[position];
    }

  }; // FlipSign<>

  template <typename T>
  class Optimize : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

    int opt_param;

    Optimize(int param) : opt_param(param) {}

  protected:
    virtual void actionOnSparseVector(SparseVector_t& v) const override
    {
      if (opt_param < 0)
        v.optimize();
      else
        v.optimize(opt_param);
    }

    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "optimize the sparse vector (";
      if (opt_param < 0)
        out << "default settings";
      else
        out << "optimization parameter: " << opt_param;
      out << ")";
    }

  }; // Optimize<>

  template <typename T>
  class FailTest : public BaseAction<T> {
  public:
    using Base_t = BaseAction<T>;
    using typename Base_t::Data_t;
    using typename Base_t::SparseVector_t;
    using typename Base_t::TestClass_t;
    using typename Base_t::Vector_t;

  protected:
    virtual void actionOnVector(Vector_t& v) const override { v.push_back(Data_t(v.size())); }

    virtual void doDescribe(TestClass_t&, std::ostream& out) const override
    {
      out << "designed failure: changes only vector";
    }

  }; // FailTest<>

} // namespace actions

//------------------------------------------------------------------------------
//---  Template implementation
//------------------------------------------------------------------------------

template <typename T>
template <typename ITER>
ITER actions::BaseAction<T>::skipRange(ITER start, ITER const end)
{
  for (; start != end; ++start)
    if (isVoid(*start)) break;
  return start;
} // actions::BaseAction::skipRange()

template <typename T>
template <typename ITER>
ITER actions::BaseAction<T>::skipVoid(ITER start, ITER const end)
{
  for (; start != end; ++start)
    if (!isVoid(*start)) break;
  return start;
} // actions::BaseAction::skipVoid()

template <typename T>
template <typename ITER>
ITER actions::BaseAction<T>::rangeStart(ITER const begin, ITER last)
{
  if (last == begin) return begin;
  while (--last != begin)
    if (isVoid(*last)) break;
  return std::next(last);
} // actions::BaseAction::rangeStart()

template <typename T>
template <typename ITER>
ITER actions::BaseAction<T>::voidStart(ITER const begin, ITER last)
{
  if (last == begin) return begin;
  while (--last != begin)
    if (!isVoid(*last)) break;
  return std::next(last);
} // actions::BaseAction::voidStart()

template <typename T>
template <typename ITER>
ITER actions::BaseAction<T>::findRangeStart(Vector_t const& v, std::size_t i, ITER start)
{
  auto const end = v.cend();

  start = skipVoid(start, end);
  while (start != end) {
    if (i-- == 0) break;
    start = skipRange(start, end);
    start = skipVoid(start, end);
  } // while
  return start;
} // actions::BaseAction::findRangeStart()

template <typename T>
template <typename ITER>
ITER actions::BaseAction<T>::findVoidStart(Vector_t const& v, std::size_t i, ITER start)
{
  auto const end = v.cend();

  start = skipRange(start, end);
  while (start != end) {
    if (i-- == 0) break;
    start = skipVoid(start, end);
    start = skipRange(start, end);
  } // while
  return start;
} // actions::BaseAction::findVoidStart()

//------------------------------------------------------------------------------

/// A simple test suite
int main()
{

  // let's test on:
  typedef float Data_t;

  //  typedef std::vector<Data_t> Vector_t;
  //  typedef lar::sparse_vector<Data_t> SparseVector_t;

  // first test: instanciate
  TestManagerClass<Data_t> Test;

  Test(actions::BaseAction<Data_t>());

  Test(actions::Resize<Data_t>(10));

  Test(actions::Assign<Data_t>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));

  Test(actions::Clear<Data_t>());

  Test(actions::AssignMove<Data_t>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));

  Test(actions::Clear<Data_t>());

  Test(actions::Insert<Data_t>(5, {5, 6, 7}));

  Test(actions::Insert<Data_t>(15, {15, 16, 17}));

  Test(actions::Insert<Data_t>(11, {11, 12, 13}));

  Test(actions::Insert<Data_t>(3, {3}));

  Data_t new_size = Data_t(Test.current_vector_size());
  for (Data_t value :
       {new_size++, new_size++ * 0, new_size++, new_size++, new_size++ * 0, new_size++})
    Test(actions::PushBack<Data_t>(value));

  Test(actions::ResizeWith<Data_t>(new_size += 2, -1));

  Test(actions::Resize<Data_t>(new_size += 5));

  Test(actions::Insert<Data_t>(7, {7, 8}));

  Test(actions::Insert<Data_t>(18, {18, 19, 20}));

  Test(actions::Insert<Data_t>(4, {4, 5}));

  Test(actions::Insert<Data_t>(22, {22}));

  Test(actions::Insert<Data_t>(19, {19}));

  Test(actions::Erase<Data_t>(7, 17));

  Test(actions::Erase<Data_t>(21, 23));

  Test(actions::PrintSparseVector<Data_t>());

  for (size_t i = 0; i < Test.current_vector_size(); ++i)
    Test.run(actions::PrintRange<Data_t>(i));

  Test(actions::Insert<Data_t>(15, {15, 16, 17, 18, 19, 20, 21, 22, 23}));

  Test(actions::Insert<Data_t>(9, {9, 10}));

  Test(actions::Erase<Data_t>(9, 11));

  Test(actions::Insert<Data_t>(9, {9, 10}));

  Test(actions::EraseRangeAt<Data_t>(11));

  Test(actions::EraseRangeAt<Data_t>(10));

  Test(actions::SetValue<Data_t>(15, -15.));

  Test(actions::FlipSign<Data_t>(15));

  for (size_t i = 13; i < 16; ++i)
    Test(actions::SetElement<Data_t>(i, i));

  Test(actions::Scale(15U, 4.0F));

  Test(actions::Scale(15U, 0.25F));

  Test(actions::ScaleAll(2.0F));

  Test(actions::ScaleAll(0.5F));

  Test(actions::UnsetElement<Data_t>(14));

  Test(actions::UnsetElement<Data_t>(15));

  Test(actions::UnsetElement<Data_t>(13));

  Test(actions::UnsetElement<Data_t>(26));

  Test(actions::UnsetElement<Data_t>(25));

  Test(actions::PrintNonVoid<Data_t>());

  Test(actions::Optimize<Data_t>(-1));

  // at this point:
  // (31) [2] {
  //      0     0    0 [  3    4    5    6 ]  0     0     0
  //      0     0    0    0    0    0 [ 16   17    18    19
  //     20    21   22   23   -1 ]  0    0    0     0     0
  //      0
  //   }

  Test(actions::Add<Data_t>(5, {7, 8, 7, 8}, 10));
  // at this point:
  // (31) [2] {
  //      0     0    0 [  3    4   12   13   17    18 ]   0
  //      0     0    0    0    0    0 [ 16   17    18    19
  //     20    21   22   23   -1 ]  0    0    0     0     0
  //      0
  //   }

  Test(actions::Add<Data_t>(5, {20, 20, 20, 20, 8, 7, 8, 7, 8, 7, 8, 7}, 30));
  // at this point:
  // (31) [2] {
  //      0     0    0 [  3    4   32   33   37    38    38
  //     37    38   37   38   37   38   23   17    18    19
  //     20    21   22   23   -1 ]  0    0    0     0     0
  //      0
  //   }

  Test(actions::Add<Data_t>(27,
                            {
                              7,
                              8,
                            },
                            20));
  // at this point:
  // (31) [2] {
  //      0     0    0 [  3    4   32   33   37    38    38
  //     37    38   37   38   37   38   23   17    18    19
  //     20    21   22   23   -1 ]  0    0 [ 27    28 ]   0
  //      0
  //   }

  Test(actions::Add<Data_t>(4, TestManagerClass<Data_t>::Vector_t(22, 10.0), 15.0));
  // at this point:
  // (31) [2] {
  //      0     0    0 [  3   14   42   43   47    48    48
  //     47    48   47   48   47   48   33   27    28    29
  //     30    31   32   33    9   25 ]  0 [ 27    28 ]   0
  //      0
  //   }

  Test(actions::Add<Data_t>(2, TestManagerClass<Data_t>::Vector_t(22, -10.0), 12.0));
  // at this point:
  // (31) [2] {
  //      0     0 [  2    3    4   32   33   37    38    38
  //     37    38   37   38   37   38   23   17    18    19
  //     20    21   22   23   -1   25 ]  0 [ 27    28 ]   0
  //      0
  //   }

  Test(actions::Add<Data_t>(26, {16, 0}, 10.0));
  // at this point:
  // (31) [2] {
  //      0     0 [  2    3    4   32   33   37    38    38
  //     37    38   37   38   37   38   23   17    18    19
  //     20    21   22   23   -1   25   26   27    28 ]   0
  //      0
  //   }

  Test(actions::Truncate<Data_t>(new_size -= 3));

  Test(actions::Truncate<Data_t>(16));

  Test(actions::Truncate<Data_t>(15));

  Test(actions::Truncate<Data_t>(14));

  Test(actions::Truncate<Data_t>(9));

  Test(actions::Truncate<Data_t>(3));

  Test(actions::Truncate<Data_t>(2));

#ifdef SPARSE_VECTOR_TEST_FAIL
  // enable this to verify that the error detection works
  Test(actions::FailTest<Data_t>());
  Test.recover();
#endif // SPARSE_VECTOR_TEST_FAIL

  return Test.summary();
} // main()
