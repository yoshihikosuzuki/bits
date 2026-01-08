def ordered_repr(*fields: str):
    """Class decorator to provide an ordered __repr__ and an __init__ method.

    The generated __init__ assigns the given ordered fields from positional
    and keyword arguments, ensures annotated attributes exist (using class
    defaults when available), and calls the class's original ``__post_init__``
    if present.

    Example: @ordered_repr("name", "seq")
    """

    def deco(cls):
        def __repr__(self):
            text = ", ".join(f"{k}={getattr(self, k)!r}" for k in fields)
            return f"{self.__class__.__name__}({text})"

        cls.__repr__ = __repr__

        # Preserve any existing __post_init__ so dataclass behaviours remain
        orig_post = getattr(cls, "__post_init__", None)

        # Only generate and assign an __init__ if the class did not define one
        # itself (i.e. avoid overwriting explicit __init__ to keep signatures).
        if "__init__" not in cls.__dict__:

            def __init__(self, *args, **kwargs):
                # Assign positional args to fields in order
                for k, v in zip(fields, args):
                    setattr(self, k, v)
                # Assign keyword args for specified fields
                for k, v in kwargs.items():
                    if k in fields:
                        setattr(self, k, v)

                # Collect annotations from class MRO, resolving forward refs when possible
                try:
                    from typing import get_type_hints

                    annotations = {}
                    for base in reversed(cls.__mro__):
                        try:
                            anns = get_type_hints(base)
                        except Exception:
                            anns = getattr(base, "__annotations__", {}) or {}
                        annotations.update(anns)
                except Exception:
                    annotations = {}
                    for base in reversed(cls.__mro__):
                        annotations.update(getattr(base, "__annotations__", {}) or {})

                # Ensure annotated attributes exist (set to class default or None)
                for name in annotations:
                    if not hasattr(self, name):
                        if hasattr(cls, name):
                            setattr(self, name, getattr(cls, name))
                        else:
                            setattr(self, name, None)

                # Call original __post_init__ if present (dataclass semantics)
                if orig_post:
                    orig_post(self)

            cls.__init__ = __init__

            # Create a nicer signature for inspection tools like IPython `?`.
            try:
                import inspect

                # Collect annotations from class MRO, resolving forward refs when possible
                try:
                    from typing import get_type_hints

                    annotations = {}
                    for base in reversed(cls.__mro__):
                        try:
                            anns = get_type_hints(base)
                        except Exception:
                            anns = getattr(base, "__annotations__", {}) or {}
                        annotations.update(anns)
                except Exception:
                    annotations = {}
                    for base in reversed(cls.__mro__):
                        annotations.update(getattr(base, "__annotations__", {}) or {})

                params = []
                for name in fields:
                    ann = annotations.get(name, inspect._empty)
                    params.append(
                        inspect.Parameter(
                            name,
                            inspect.Parameter.POSITIONAL_OR_KEYWORD,
                            annotation=ann,
                        )
                    )
                sig = inspect.Signature(parameters=params)
                cls.__init__.__signature__ = sig
                cls.__signature__ = sig

                # Populate __annotations__ on the generated __init__ so help() shows types
                init_ann = {
                    name: annotations.get(name, inspect._empty) for name in fields
                }
                init_ann["return"] = None
                try:
                    cls.__init__.__annotations__ = init_ann
                except Exception:
                    pass
            except Exception:
                pass
        return cls

    return deco


def verbose_repr(*fixed_fields: str):
    """Class decorator to provide a verbose __repr__ that lists fixed fields
    first and then any additional attributes present on the instance.

    Example: @verbose_repr("chrom", "b", "e")
    """

    def deco(cls):
        def __repr__(self):
            fixed = list(fixed_fields)
            fixed_vals = [f"{name}={repr(getattr(self, name))}" for name in fixed]
            other = list(set(vars(self).keys()) - set(fixed))
            other_vals = [f"{name}={repr(getattr(self, name))}" for name in other]
            text = ", ".join(fixed_vals + other_vals)
            return f"{self.__class__.__name__}({text})"

        cls.__repr__ = __repr__
        return cls

    return deco
