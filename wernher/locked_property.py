import inspect
import numpy as np

class LockError(Exception):
    def __init__(self,call_stack=None):
        if call_stack is not None:
            self.call_stack = call_stack[:]
        else:
            self.call_stack = None

    def __str__(self):
        msg  = 'Not enough information.'
        if self.call_stack is not None:
            msg += '\n  call stack:\n    '
            msg += '\n    '.join(self.call_stack)
        return msg

class OverspecifiedError(Exception):
    def __init__(self,property_name):
        self.property_name = property_name

    def __str__(self):
        return 'Overspecified with \''+self.property_name+'\'.'

class cached_property(object):
    """
    A read-only @property that is only evaluated once.
    """
    def __init__(self, fget, doc=None):
        self.__doc__ = doc or fget.__doc__
        self.__name__ = fget.__name__
        self.fget = fget

    def __get__(self, obj, cls):
        if obj is None:
            return self
        if not hasattr(obj,'_calculated_properties'):
            obj._calculated_properties = list()
        if self.__name__ not in obj.__dict__:
            obj.__dict__[self.__name__] = self.fget(obj)
            obj._calculated_properties.append(self.__name__)
        return obj.__dict__[self.__name__]

    def __delete__(self, obj):
        if hasattr(obj,'_calculated_properties'):
            for p in obj._calculated_properties:
                if p is not self.__name__:
                    del obj.__dict__[p]
            obj._calculated_properties.clear()
        del obj.__dict__[self.__name__]

class locked_property(cached_property):
    """
    A read-only @property that is only evaluated once.
    And can only be entered once (locked).
    """
    def __init__(self, fget, doc=None):
        self.__doc__ = doc or fget.__doc__
        self.__name__ = fget.__name__
        self.fget = self.lock(fget)

    def __set__(self, obj, val):
        if self.__name__ in obj.__dict__:
            self.__delete__(obj)
        try:
            getattr(obj,self.__name__)
            raise OverspecifiedError(self.__name__)
        except LockError:
            if hasattr(val,'__iter__'):
                val = np.asarray(val, dtype=np.float64)
            try:
                obj.__dict__[self.__name__] = np.float64(val)
            except TypeError:
                obj.__dict__[self.__name__] = val

    def lock(self,fn):
        def locked_fget(obj):
            if not hasattr(obj,'call_stack'):
                obj.call_stack = list()
            if self.__name__ in obj.call_stack:
                raise LockError(obj.call_stack)
            else:
                obj.call_stack.append(self.__name__)
                try:
                    val = fn(obj)
                    self.debug = False
                    if self.debug:
                        print(fn.__name__)
                        frame = inspect.stack()[2][0]
                        try:
                            ctx = inspect.getframeinfo(frame)
                            print('   ',ctx)
                        finally:
                            del frame

                finally:
                    obj.call_stack.remove(self.__name__)
                return val
        return locked_fget

def property_alias(name):
    def getter(self):
        return getattr(self,name)
    def setter(self,val):
        setattr(self,name,val)
    def deleter(self):
        delattr(self,name)
    return property(getter,setter,deleter)
