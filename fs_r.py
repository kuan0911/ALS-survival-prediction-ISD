import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import warnings
import numbers
from operator import attrgetter

from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis

from sklearn.model_selection import GridSearchCV, KFold
from sklearn.pipeline import make_pipeline
from sklearn.utils import safe_sqr


def coxen_fs_r(
        data,
        alpha_min_ratio=0.01,
        verbose: bool = True
):
    l1_ratio_list = [0.9, 0.7, 0.5]

    features = data.columns.tolist()
    features.remove('time')
    features.remove('event')

    X = data[features]
    event_times = data['time'].values
    event_indicators = data['event'].values.astype('bool')

    y = np.empty(dtype=[('cens', bool), ('time', np.float64)], shape=event_times.shape[0])
    y['cens'] = event_indicators
    y['time'] = event_times

    best_cindex = 0
    for l1_ratio in tqdm(l1_ratio_list, total=len(l1_ratio_list), disable=verbose):
        coxnet_pipe = make_pipeline(
            CoxnetSurvivalAnalysis(l1_ratio=l1_ratio, alpha_min_ratio=alpha_min_ratio, max_iter=100)
        )
        warnings.simplefilter("ignore", UserWarning)
        coxnet_pipe.fit(X, y)

        estimated_alphas = coxnet_pipe.named_steps["coxnetsurvivalanalysis"].alphas_
        cv = KFold(n_splits=5, shuffle=True, random_state=0)
        gcv = GridSearchCV(
            make_pipeline(CoxnetSurvivalAnalysis(l1_ratio=l1_ratio)),
            param_grid={"coxnetsurvivalanalysis__alphas": [[v] for v in estimated_alphas]},
            cv=cv,
            error_score=0.5,
            n_jobs=1).fit(X, y)

        if gcv.best_score_ > best_cindex:
            best_cindex = gcv.best_score_
            best_model = gcv.best_estimator_.named_steps["coxnetsurvivalanalysis"]
    # best_model = gcv.best_estimator_.named_steps["coxnetsurvivalanalysis"]
    # best_coefs = pd.DataFrame(
    #     best_model.coef_,
    #     index=data.columns.drop(['time', 'event']),
    #     columns=["coefficient"]
    # )

    non_zero_idx = best_model.coef_.squeeze() != 0
    features_selected = data.columns.drop(['time', 'event'])[non_zero_idx]
    fs_columns = features_selected.tolist()
    fs_columns.append('time')
    fs_columns.append('event')
    return data[fs_columns]


def rfe_fs_r(df, n_features_to_select=None, step=1, verbose=True):
    # Parameter step_score controls the calculation of self.scores_
    # step_score is not exposed to users
    # and is used when implementing RFECV
    # self.scores_ will not be calculated when calling _fit through fit
    n_features = df.shape[1] - 2
    n_features_to_select = int(n_features_to_select)

    if 0.0 < step < 1.0:
        step = int(max(1, step * n_features))
    else:
        step = int(step)
    if step <= 0:
        raise ValueError("Step must be >0")

    X = df.drop(columns=['event', 'time'], inplace=False).values
    T = df['time'].values
    E = df['event'].values.astype('bool')

    y = np.empty(dtype=[('cens', bool), ('time', np.float64)], shape=T.shape[0])
    y['cens'] = E
    y['time'] = T

    support_ = np.ones(n_features, dtype=bool)
    ranking_ = np.ones(n_features, dtype=int)

    # Elimination
    while np.sum(support_) > n_features_to_select:
        # Remaining features
        features = np.arange(n_features)[support_]

        # Rank the remaining features
        estimator = CoxPHSurvivalAnalysis()
        if verbose:
            print("Fitting estimator with %d features." % np.sum(support_))

        estimator.fit(X[:, features], y)

        # Get importance and rank them
        importances = _get_feature_importances(
            estimator,
            "auto",
            transform_func="square",
        )
        ranks = np.argsort(importances)

        # for sparse case ranks is matrix
        ranks = np.ravel(ranks)

        # Eliminate the worse features
        threshold = min(step, np.sum(support_) - n_features_to_select)

        # Compute step score on the previous selection iteration
        # because 'estimator' must use features
        # that have not been eliminated yet
        support_[features[ranks][:threshold]] = False
        ranking_[np.logical_not(support_)] += 1

    # Set final attributes
    features = np.arange(n_features)[support_]
    feature_names = df.columns.drop(['event', 'time']).to_numpy()[features]
    print("Selected features by RFE are: {}".format(feature_names))
    feature_names = feature_names.tolist()
    feature_names.append('time')
    feature_names.append('event')
    # estimator_ = CoxPHSurvivalAnalysis()
    # estimator_.fit(X[:, features], y)
    #
    # # Compute step score when only n_features_to_select features left
    # self.n_features_ = support_.sum()
    # self.support_ = support_
    # self.ranking_ = ranking_

    return df[feature_names]


def _get_feature_importances(estimator, getter, transform_func=None, norm_order=1):
    """
    Retrieve and aggregate (ndim > 1)  the feature importances
    from an estimator. Also optionally applies transformation.
    Parameters
    ----------
    estimator : estimator
        A scikit-learn-type estimator from which we want to get the feature
        importances.
    getter : "auto", str or callable
        An attribute or a callable to get the feature importance. If `"auto"`,
        `estimator` is expected to expose `coef_` or `feature_importances`.
    transform_func : {"norm", "square"}, default=None
        The transform to apply to the feature importances. By default (`None`)
        no transformation is applied.
    norm_order : int, default=1
        The norm order to apply when `transform_func="norm"`. Only applied
        when `importances.ndim > 1`.
    Returns
    -------
    importances : ndarray of shape (n_features,)
        The features importances, optionally transformed.
    """
    if isinstance(getter, str):
        if getter == "auto":
            if hasattr(estimator, "coef_"):
                getter = attrgetter("coef_")
            elif hasattr(estimator, "feature_importances_"):
                getter = attrgetter("feature_importances_")
            else:
                raise ValueError(
                    "when `importance_getter=='auto'`, the underlying "
                    f"estimator {estimator.__class__.__name__} should have "
                    "`coef_` or `feature_importances_` attribute. Either "
                    "pass a fitted estimator to feature selector or call fit "
                    "before calling transform."
                )
        else:
            getter = attrgetter(getter)
    elif not callable(getter):
        raise ValueError("`importance_getter` has to be a string or `callable`")
    importances = getter(estimator)

    if transform_func is None:
        return importances
    elif transform_func == "norm":
        if importances.ndim == 1:
            importances = np.abs(importances)
        else:
            importances = np.linalg.norm(importances, axis=0, ord=norm_order)
    elif transform_func == "square":
        if importances.ndim == 1:
            importances = safe_sqr(importances)
        else:
            importances = safe_sqr(importances).sum(axis=0)
    else:
        raise ValueError(
            "Valid values for `transform_func` are "
            + "None, 'norm' and 'square'. Those two "
            + "transformation are only supported now"
        )

    return importances


