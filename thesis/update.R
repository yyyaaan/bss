p14, l-5: We usually say that the data are whitened, not matrices. So first whiten the data using R_0, then computed R_\tau using whitened data.

about whitenend data vs whitened matrices: it seems that I cannot use whitented data for JADE. This is becasue the whitening applies only to the decomposed autocovariance, and thus is not achievable by whitening data directly.

f