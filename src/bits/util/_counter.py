from collections import Counter


class RelCounter(Counter):
    """Subclass of Counter with a method returning the relative frequencies."""

    def relative(self):
        tot = sum(self.values())
        return {k: v / tot * 100 for k, v in self.items()}
