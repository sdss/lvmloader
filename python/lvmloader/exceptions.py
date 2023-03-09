# !usr/bin/env python
# -*- coding: utf-8 -*-
#
# Licensed under a 3-clause BSD license.
#
# @Author: Brian Cherinka
# @Date:   2017-12-05 12:01:21
# @Last modified by:   Brian Cherinka
# @Last Modified time: 2017-12-05 12:19:32

from __future__ import absolute_import, division, print_function


class LvmloaderError(Exception):
    """A custom core Lvmloader exception"""

    def __init__(self, message=None):

        message = 'There has been an error' \
            if not message else message

        super(LvmloaderError, self).__init__(message)


class LvmloaderNotImplemented(LvmloaderError):
    """A custom exception for not yet implemented features."""

    def __init__(self, message=None):

        message = 'This feature is not implemented yet.' \
            if not message else message

        super(LvmloaderNotImplemented, self).__init__(message)


class LvmloaderAPIError(LvmloaderError):
    """A custom exception for API errors"""

    def __init__(self, message=None):
        if not message:
            message = 'Error with Http Response from Lvmloader API'
        else:
            message = 'Http response error from Lvmloader API. {0}'.format(message)

        super(LvmloaderAPIError, self).__init__(message)


class LvmloaderApiAuthError(LvmloaderAPIError):
    """A custom exception for API authentication errors"""
    pass


class LvmloaderMissingDependency(LvmloaderError):
    """A custom exception for missing dependencies."""
    pass


class LvmloaderWarning(Warning):
    """Base warning for Lvmloader."""


class LvmloaderUserWarning(UserWarning, LvmloaderWarning):
    """The primary warning class."""
    pass


class LvmloaderSkippedTestWarning(LvmloaderUserWarning):
    """A warning for when a test is skipped."""
    pass


class LvmloaderDeprecationWarning(LvmloaderUserWarning):
    """A warning for deprecated features."""
    pass
