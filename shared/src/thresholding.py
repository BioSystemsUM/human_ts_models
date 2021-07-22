from numpy import log


def global_thresholding(sample_series, gtu, gtl, lt, maxexp):
    return (sample_series / gtu).apply(log).clip(-maxexp, maxexp).to_dict()


def local1_thresholding(sample_series, gtu, gtl, lt, maxexp):
    activity = (sample_series / gtu).apply(log).clip(-maxexp, maxexp)
    gt_active = activity >= 0
    activity[gt_active] = (log(sample_series[gt_active] / lt[gt_active]) * maxexp).clip(-maxexp, maxexp)
    return activity.to_dict()


def local2_thresholding(sample_series, gtu, gtl, lt, maxexp):
    upp_activity = (1 + (sample_series / gtu).apply(log)).clip(-maxexp, 1 + maxexp)

    gtu_inactive = upp_activity < 1
    low_activity = (sample_series / gtl).apply(log).clip(-maxexp, maxexp)
    gtl_maybes, gtl_lows = (low_activity >= 0) & gtu_inactive, low_activity < 0
    upp_activity[gtl_lows] = low_activity[gtl_lows]
    activity_maybe = (sample_series[gtl_maybes] / lt[gtl_maybes]). \
        apply(log).clip(-maxexp, maxexp)
    upp_activity[gtl_maybes] = activity_maybe.clip(-1, 1)
    return upp_activity.to_dict()

THRESHOLDING_FUNCTIONS = {
    'global': global_thresholding,
    'local1': local1_thresholding,
    'local2': local2_thresholding
}