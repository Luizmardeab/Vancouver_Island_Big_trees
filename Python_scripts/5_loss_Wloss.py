import tensorflow as tf
import tensorflow.keras.backend as K
import tensorflow_probability as tfp


# Precomputed bin edges and bin weights for Height
bin_edges = tf.constant([i for i in range(96)], dtype=tf.float32)  # Bins from 0 to 96
bin_weights = tf.constant([0.0081314 , 0.01731114, 0.01806327, 0.01952283, 0.01989786,
       0.01999887, 0.01984188, 0.01993265, 0.02026657, 0.02022913,
       0.01984085, 0.01915748, 0.01907826, 0.01871557, 0.01830211,
       0.01791165, 0.01757751, 0.01725701, 0.01697489, 0.01673199,
       0.01650983, 0.01632886, 0.01616751, 0.01604569, 0.01595099,
       0.0159067 , 0.01584951, 0.01584459, 0.01585675, 0.01591508,
       0.01600977, 0.01613661, 0.01633729, 0.01656966, 0.0168444 ,
       0.01718832, 0.01761204, 0.01812263, 0.01872677, 0.01942361,
       0.02020101, 0.02104552, 0.0220224 , 0.02306922, 0.02427871,
       0.02563085, 0.02715977, 0.02894054, 0.03093484, 0.03330584,
       0.0360381 , 0.0393531 , 0.04317577, 0.04781715, 0.05313619,
       0.05940371, 0.06658575, 0.07503753, 0.08474138, 0.09529222,
       0.10795234, 0.12173416, 0.13795952, 0.1544682 , 0.17299899,
       0.193213  , 0.21440476, 0.23658717, 0.26098365, 0.28788931,
       0.31377732, 0.34064878, 0.37021805, 0.40956099, 0.43725275,
       0.47415203, 0.51505353, 0.5562427 , 0.60350569, 0.63397891,
       0.67055223, 0.71312238, 0.7302578 , 0.79818857, 0.81890223,
       0.85348593, 0.87678828, 0.89172047, 0.88946008, 0.9552009 ,
       0.98513801, 0.95589889, 1.00240867, 1.00727871, 0.97165259], dtype=tf.float32)

def assign_bin_weights(y_true):
    """Assign precomputed bin weights based on canopy height bins."""
    mask_gte_zero = K.greater_equal(y_true, 0) #Mask values bellow 0
    mask_lt_percentile_95 = K.less(y_true, 96) #Mask values above 96
    mask = K.all([mask_gte_zero, mask_lt_percentile_95], axis=0)  #Create mask
    y_true_masked = tf.boolean_mask(y_true, mask)
    y_true_flat = tf.reshape(y_true_masked, [-1])  # Flatten
    bin_indices = tfp.stats.find_bins(y_true_flat, bin_edges)  # Find bins
    bin_indices = tf.cast(bin_indices, tf.int32)  # Convert to int32
    sample_weights = tf.gather(bin_weights, bin_indices)  # Assign weights
    return tf.reshape(sample_weights, tf.shape(y_true_masked))  # Reshape to match input shape

def weighted_loss(y_true, y_pred, loss_fn):
    """Apply precomputed bin-based weighting to a loss function."""
    weights = assign_bin_weights(y_true)  # Get weights for this batch
    mask_gte_zero = K.greater_equal(y_true, 0) #Mask values bellow 0
    mask_lt_percentile_95 = K.less(y_true, 96) #Mask values above 96
    mask = K.all([mask_gte_zero, mask_lt_percentile_95], axis=0)  #Create mask
    y_true_masked = tf.boolean_mask(y_true, mask)
    y_pred_masked = tf.boolean_mask(y_pred, mask)
    loss = loss_fn(y_true_masked, y_pred_masked)  # Compute base loss
    return tf.reduce_mean(loss * weights)  # Apply weighting

def gaussian_nll_loss(y_true, y_pred):
    """Gaussian NLL loss with bin-based weights and valid-value masking."""
    # Ensure y_true has shape (B, H, W)
    if y_true.shape.rank == 4:
        y_true = tf.squeeze(y_true, axis=-1)  # From (B, H, W, 1) to (B, H, W)
    weights = assign_bin_weights(y_true)  # Shape: (B, H, W)
    
    # Create mask to exclude invalid values
    mask = tf.logical_and(
        tf.greater_equal(y_true, 0.0),
        tf.less_equal(y_true, 96.0)
    )  # Shape: (B, H, W)

    # Extract mean and log_var from y_pred
    mean = y_pred[..., 0]       # Shape: (B, H, W)
    log_var = y_pred[..., 1]    # Shape: (B, H, W)

    # Apply mask to all tensors
    y_true_masked = tf.boolean_mask(y_true, mask)
    mean_masked = tf.boolean_mask(mean, mask)
    log_var_masked = tf.boolean_mask(log_var, mask)
    # weights = tf.boolean_mask(weights, mask)

    # Clip log_var for numerical stability
    log_var_masked = tf.clip_by_value(log_var_masked, -10.0, 10.0)
    precision = tf.exp(-log_var_masked)

    # Gaussian NLL formula
    loss = 0.5 * (tf.math.log(2.0 * np.pi) + log_var_masked +
                  tf.square(y_true_masked - mean_masked) * precision)

    return tf.reduce_mean(loss * weights)

def laplace_loss(y_true, y_pred):
    """Modified Laplacian loss with sample weighting."""
    return weighted_loss(y_true, y_pred, lambda y_t, y_p: tf.abs(y_t - y_p))

def rmse(y_true, y_pred):
    """Modified RMSE with sample weighting."""
    return K.sqrt(weighted_loss(y_true, y_pred, lambda y_t, y_p: K.square(y_t - y_p)))

def normalized_rmse(y_true, y_pred):
    """Modified Normalized RMSE with sample weighting."""
    weights = assign_bin_weights(y_true) #Calculate weights
    mask_gte_zero = K.greater_equal(y_true, 0) #Mask values bellow 0
    mask_lt_percentile_95 = K.less(y_true, 96) #Mask values above 96
    mask = K.all([mask_gte_zero, mask_lt_percentile_95], axis=0)  #Create mask
    y_true_masked = tf.boolean_mask(y_true, mask)*weights
    y_pred_masked = tf.boolean_mask(y_pred, mask)*weights
    rmse_value = rmse(y_true_masked, y_pred_masked)
    mean_true = K.mean(y_true_masked)
    return (rmse_value / (mean_true + K.epsilon())) * 100

def mae(y_true, y_pred):
    return K.mean(K.abs(y_pred - y_true))


def bias(y_true, y_pred):
    """Modified Bias with sample weighting."""
    return weighted_loss(y_true, y_pred, lambda y_t, y_p: y_p - y_t)

def r_squared(y_true, y_pred):
    """Modified R-squared with sample weighting."""
    weights = assign_bin_weights(y_true)
    mask_gte_zero = K.greater_equal(y_true, 0) #Mask values bellow 0
    mask_lt_percentile_95 = K.less(y_true, 96) #Mask values above 96
    mask = K.all([mask_gte_zero, mask_lt_percentile_95], axis=0)  #Create mask
    y_true_masked = tf.boolean_mask(y_true, mask)*weights
    y_pred_masked = tf.boolean_mask(y_pred, mask)*weights

    residual = K.sum(K.square(y_true_masked - y_pred_masked))
    total = K.sum(K.square(y_true_masked - K.mean(y_true_masked)))
    
    return 1 - (residual / (total + K.epsilon()))
