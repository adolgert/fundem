Lifetable
=========

.. index:: survival, px

.. function:: first_moment_survival(mx, ax, nx)

    :param array[pop,age] mx: Mortality rate :math:`m_x`.
    :param array[pop,age] ax: A starting set of values for :math:`a_x`.
                             These could be generated with
                             ``ax = cm_mean_age(mx, nx)``. *This will
                             be changed in place.*
    :param array[age] nx: Interval sizes which are uniform for all age
                             groups.
    :return: Survival :math:`{}_np_x`.
    :rtype: array[pop,age]

    This is :math:`{}_np_x`. It is the survival for intervals.
    Given :math:`m_x`, and :math:`a_x`, this can make an exact calculation.

    .. math::

       {}_np_x = \frac{1-m_x a_x}{1+ m_x(n_x - a_x)}


.. index:: deaths, population, dx, lx

.. function:: first_moment_population(mx, ax, nx)

    :param array[pop,age] mx: Mortality rate :math:`m_x`.
    :param array[pop,age] ax: A starting set of values for :math:`a_x`.
                             These could be generated with
                             ``ax = cm_mean_age(mx, nx)``. *This will
                             be changed in place.*
    :param array[age] nx: Interval sizes which are uniform for all age
                             groups.
    :return: Survival :math:`(l_x, {}_nd_x)`.
    :rtype: (array[pop,age],array[pop,age])

    Calculates both :math:`(l_x, {}_nd_x)` from mortality rate and mean
    age. This is an exact calculation.

    .. math::

       {}_np_x = \frac{1-m_x a_x}{1+ m_x(n_x - a_x)}


.. index:: life expectancy, LE

.. function:: first_moment_period_life_expectancy(mx, ax, nx)

    :param array[pop,age] mx: Mortality rate :math:`m_x`.
    :param array[pop,age] ax: A starting set of values for :math:`a_x`.
                             These could be generated with
                             ``ax = cm_mean_age(mx, nx)``. *This will
                             be changed in place.*
    :param array[age] nx: Interval sizes which are uniform for all age
                             groups.
    :return: Period life expectancy :math:`\mathring{e}_x`.
    :rtype: array[pop,age]

    Calculates life expectancy using recursive method.
    This uses the mean age of death to specify the first moment of the
    distribution of deaths. The period life expectancy is a statistic
    that uses the mortality rate for all populations in the same year.
    This is complete, not curtate.

    This uses a recursive algorithm,

    .. math::

       \mathring{e}_x=\frac{n_x + (1-m_x\:{}_na_x)\mathring{e}_{x+n_x}}
       {1+m_x(n_x-{}_na_x)}

    starting from :math:`\mathring{e}_x=a_x` for the last, half-open interval.


.. index:: graduation method

.. function:: graduation_method(mx, nx, ax)

    :param array[pop,age] mx: Mortality rate :math:`m_x`.
    :param array[age] nx: Interval sizes which are uniform for all age
                             groups.
    :param array[pop,age] ax: A starting set of values for :math:`a_x`.
                             These could be generated with
                             ``ax = cm_mean_age(mx, nx)``. *This will
                             be changed in place.*

    This is an iterative method to determine :math:`a_x` from :math:`m_x`,
    which Preston calls the *graduation method.* It estimates an initial
    value for :math:`a_x` and then smooths it using a cubic spline.
    The equation in Preston is

    .. math::

       {}_na_x = \frac{-\frac{n}{24}{}_nd_{x-n}+\frac{n}{2}{}_nd_x +
       \frac{n}{24}{}_nd_{x+n}}{{}_nd_x}.

    The intervals :math:`n_x` must be equal (usually 5 years), and this
    can only smooth :math:`a_x` for intervals which have left
    and right death counts. Preston defines this in terms of deaths,
    but we implement it against mortality, just multiplying by :math:`l_x`.
