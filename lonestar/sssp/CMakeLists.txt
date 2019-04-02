app(sssp SSSP.cpp bf/bloom_filter/basic.cpp bf/hash.cpp bf/bitvector.cpp)

add_test_scale(small sssp "${BASEINPUT}/structured/rome99.gr")
add_test_scale(web sssp "${BASEINPUT}/random/r4-2e26.gr" -delta 8)
