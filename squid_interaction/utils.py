import tensorflow as tf


class SquidWrapper(tf.keras.Model):
    """
    A small wrapper class to let CLIPNET predict on one-hot encoded sequences.
    """

    def __init__(self, base_model):
        super(SquidWrapper, self).__init__()
        self.base_model = base_model

    def call(self, inputs, training=False):
        # Multiply inputs by 2 before passing to the base model
        return self.base_model(inputs * 2, training=training)


class QuantityEnsembler(tf.keras.Model):
    """
    A small wrapper class to average a bunch of models that have been passed through SquidWrapper.
    """

    def __init__(self, a_list_of_models):
        super(QuantityEnsembler, self).__init__()
        self.base_models = a_list_of_models

    def call(self, inputs, training=False):
        return tf.expand_dims(
            tf.keras.layers.Average()([m(inputs)[1] for m in self.base_models]), axis=-1
        )


class ProfileEnsembler(tf.keras.Model):
    """
    A small wrapper class to average a bunch of models that have been passed through SquidWrapper.
    """

    def __init__(self, a_list_of_models):
        super(ProfileEnsembler, self).__init__()
        self.base_models = a_list_of_models

    def call(self, inputs, training=False):
        return tf.expand_dims(
            tf.keras.layers.Average()([m(inputs)[0] for m in self.base_models]), axis=-1
        )


def corr(x, y, pseudocount=1e-6):
    """
    Custom objective function that needs to be included when loading
    CLIPNET instances.
    """
    mx = tf.math.reduce_mean(x)
    my = tf.math.reduce_mean(y)
    xm, ym = x - mx, y - my
    num = tf.math.reduce_mean(tf.multiply(xm, ym))
    den = tf.math.reduce_std(xm) * tf.math.reduce_std(ym) + pseudocount
    r = tf.math.maximum(tf.math.minimum(num / den, 1), -1)
    return r


def softmax_sum(x, save_dir=None):
    softmax = tf.keras.layers.Softmax()
    return tf.reduce_sum(tf.stop_gradient(softmax(x)) * x, axis=-1)
