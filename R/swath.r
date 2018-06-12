#' Create a sequence of swath intervals for a given axis and increment size.
#'
#' Swath statistics are generated for slices through the spatial data.
#' \code{swath_seq} generates slice intervals based on some input spatial
#' data, a block size, and slice width.
#''
#' @author Alex Trueman
#' @param axis numeric vector of single axis coordinates (e.g., model$x).
#' @param slice_size numeric swath slice size for the axis.
#' @param block_inc numeric size of parent blocks in the model.
#'    The minimum and maximum coordinate of slicing is rounded to the
#'    nearest \code{block_inc}.
#'
#' @return numeric vector with a sequence of slice coordinates.
#' @export
#' @importFrom plyr round_any
#'
swath_seq <- function(axis, slice_size, block_inc) {

  # Calculate number of slices based on range of the axis.
  sl_num <- round(round_any(max(axis) - min(axis), block_inc) / slice_size) + 1

  # Calculate minimum and maximum slice.
  sl_min <- round_any(min(axis) - slice_size, block_inc)
  sl_max <- sl_min + (slice_size * sl_num)

  # Return numeric sequence of intervals.
  return(as.numeric(seq(sl_min, sl_max, slice_size)))

}

#' Extract the midpoint of an interval assigned to a data frame.
#'
#' The mid-point is used as an attractive label for plot axes.
#'
#' @author Alex Trueman
#' @param x is a data-frame column with intervals typically applied using the
#'   base-r `cut` function.
#' @param dp is the number of decimal places for the base-r `round` function.
#'
#' @return an atomic vector of type double.
interval_mid <- function(x, dp = 1) {

  # Extract lower and upper bounds of the intervals.
  lower <- as.double(gsub(",.*", "", gsub("\\(|\\[|\\)|\\]", "", x)))
  upper <- as.double(gsub(".*,", "", gsub("\\(|\\[|\\)|\\]", "", x)))

  # Return midpoint.
  return(as.double(round(lower + (upper - lower) / 2, dp)))

}

#' Generate means by geographical slices (swaths).
#'
#' For a given set of coordinate intervals calculate the mean of the data
#' within each interval (swath).
#' @author Alex M Trueman
#'
#' @param df Data frame with at least one coordinate field.
#' @param value Grade column being evaluated.
#' @param group Grouping column (e.g., domain).
#' @param axis Axis column for slicing (e.g., x).
#' @param slices List of slice interval coordinates as numeric vector.
#'    Can be generated from function \code{swath_seq}.
#'
#' @return data frame
#' @export
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang enquo
swath_data <- function(df, value, group, axis, slices) {

    value <- enquo(value)
    group <- enquo(group)
    axis <- enquo(axis)

    # Generate swath data by group and slice.
    df %>%
        group_by(!! group, slice = interval_mid(cut(!! axis, slices), 1)) %>%
        filter(!is.na(slice)) %>%
        summarise(n = n(), statistic = mean(!! value, na.rm = TRUE)) %>%
        select(!! group, slice, n, statistic) %>%
        arrange(!! group, slice)

}

#' Generate means by geographical slices (swaths).
#' Extended to include lower and upper confidence limits by bootstrapping.
#'
#' For a given set of coordinate intervals calculate the mean of the data
#' within each interval (swath).
#'
#' @author Alex M Trueman
#' @param df Data frame with at least one coordinate field.
#' @param value Column name to be evaluated (e.g., au).
#' @param group Grouping column (e.g., domain).
#' @param axis Axis column for slicing (e.g., x).
#' @param slices List of slice starting coordinates as numeric vector.
#'    Can be generated from function \code{swath_seq}.
#' @param reps (default = 10000) Number of bootstrap repetitions.
#' @param conf (default = 0.95) Confidence level.
#'
#' @return dataframe
#' @export
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang enquo quo_name UQ
#' @importFrom mstats boot_mean_ci
swath_data_ci <- function(df, value, group, axis, slices, reps = 10000, conf = 0.95) {

    # Set up passed arguments.
    value <- enquo(value)
    group <- enquo(group)
    group_str <- quo_name(group)
    axis <- enquo(axis)

    # Create slice field and remove records that will make bootstrap fail.
    d <- df %>%
        mutate(slice = interval_mid(cut(!! axis, slices), 1)) %>%
        group_by(!! group, slice) %>%
        mutate(n = n(), isunique = length(unique(!! value))) %>% # Count unique values per group.
        ungroup() %>%
        filter(n > 1, isunique > 1) # Make sure each group has more than 1 unique value
    # otherwise the boot.ci will fail with error.

    # Bootstrap for confidence intervals.
    mci <- boot_mean_ci(d, UQ(value), reps = reps, conf = conf, UQ(group), slice)

    # Join dataframes for plotting.
    fin <- d %>%
        group_by(!! group, slice) %>%
        summarise(n = n()) %>%
        left_join(mci, by = c(group_str, "slice")) %>%
        arrange(!! group, slice)

    return(fin)

}
