# test raw data compression
cet_test(raw_test USE_BOOST_UNIT
  LIBRARIES PRIVATE
  lardataobj::RawData
  larcoreobj::headers
)

# test data products
cet_test(RawDigit_test USE_BOOST_UNIT
  LIBRARIES PRIVATE
  lardataobj::RawData
  larcoreobj::headers
)

install_headers()
install_source()
