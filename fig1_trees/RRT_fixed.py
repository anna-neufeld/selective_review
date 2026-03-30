def run_rrt_node_inference(Xpyth, ypyth, sdy, seed, noise_sd=1.0, alpha=0.1, maxdepth=2):
  
    import numpy as np
    import pandas as pd
    from sklearn.metrics import adjusted_rand_score
    import sys
    import traceback

    from CART import RegressionTree

    np.random.seed(seed)

    # Convert R inputs to numpy
    Xpyth = np.asarray(Xpyth)
    ypyth = np.asarray(ypyth).reshape(-1)
    #sd = float(sd)
    #noise_sd = float(noise_sd)

    # ---- fit regression tree ----
    reg_tree = RegressionTree(
        max_depth=maxdepth,
        min_proportion=0,
        min_samples_split=2,
        min_bucket=1
    )

    ### I am passing in a noise_sd that has already been scaled by sd_y
    reg_tree.fit(Xpyth, ypyth, sd=noise_sd)

    muhat = reg_tree.predict(Xpyth)

    node_results = []
    nnodes = len(reg_tree.terminal_nodes)
    for node_idx in range(nnodes):
        try:
            (
                pval, dist, contrast, norm_contrast, obs_tar,
                logW, suff, sel_probs, ref_hat_layer
            ) = reg_tree.node_inference(
                node=reg_tree.terminal_nodes[node_idx],
                ngrid=10000,
                ncoarse=50,
                grid_w_const=10,
                #reduced_dim=None,
                #reduced_dim=None,
                sd=sdy, ## DO I PUT NOISE SD IN HERE TOO?
                use_cvxpy=True
            )
 
            
            #norm_contrast = np.asarray(norm_contrast, dtype=float).reshape(-1)
            norm_contrast = np.asarray(norm_contrast).reshape(-1)
            ypyth = np.asarray(ypyth).reshape(-1)
            selective_CI =  dist.equal_tailed_interval(
                observed=float(np.dot(norm_contrast, ypyth)),
                alpha=alpha
            )
            selective_CI = np.array(selective_CI)
            selective_CI *= np.linalg.norm(contrast) * sdy
            success = True

        except Exception:
            #traceback.print_exc()

            # NA-like outputs on failure
            n = Xpyth.shape[0]
            contrast = np.full(shape=(n,), fill_value=np.nan, dtype=float)
            selective_CI = np.nan, np.nan
            success = False

        ## It seems to me that this is correct
        node_results.append({
            "node_index": int(node_idx),
            "lower":  selective_CI[0],
            "upper": selective_CI[1],
            "norm_contrast": contrast, 
            "success": success
        })

    return {
        "muhat": muhat, #np.asarray(muhat, dtype=float).reshape(-1),
        "node_results": node_results
    }
