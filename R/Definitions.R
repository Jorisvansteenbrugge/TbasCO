#' Defines the S3 class object
#'
setClass('General_features', representation(name                = 'character',
                                            bins                = 'vector',
                                            sample_columns      = 'vector',
                                            rank_columns        = 'vector',
                                            KO_presence_absence = 'matrix'),
         prototype(name = 'matrix_features'))
