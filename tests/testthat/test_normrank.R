test_that("Test Create.Rank.Columns()", {

            x <- "Sample1;Bin\n10;5\n30;5\n5;5\n19;5\n90;4\n100;4"
            data  <- read.csv2(text = gsub("(?<=[a-z])\\s+", "\n",
                                           x, perl = TRUE),
                               header = T)
            features <- list('bins' = c(4, 5),
                             'sample.columns' = c(1),
                             'rank.columns'   = c(3)
                             )

            expect_equal(Create.Rank.Columns(data, features)[, 3],
                         c(0.75 , 0.25 , 1.00 , 0.50 , 1.00 , 0.50))
          })
