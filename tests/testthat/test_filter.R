test_that("Test Filtering(MAD)",
          {x <- "Gene;Sample1;Sample2;Sample3;Sample4;Sample5;Sample6;Annotation;Bin\nCAP2UW1_00001;71;36;113;36;229;101;;39\nCAP2UW1_00001;0;0;0;0;0;0;;39"
          data  <- read.csv2(text=gsub("(?<=[a-z])\\s+", "\n", x, perl=TRUE), header=T)
          expect_equal(nrow(Filter(data, sample.columns = 2:7,
                                   filter.method = 'MAD')),
                       1)})


