# sparse_vector_test tests pure header libraries
cet_test(sparse_vector_test)

# LazyVector_test tests pure header libraries
cet_test(LazyVector_test USE_BOOST_UNIT)

# flagset_test tests pure header libraries
cet_test(FlagSet_test USE_BOOST_UNIT)

install_source()
