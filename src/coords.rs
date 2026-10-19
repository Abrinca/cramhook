//! Distinct coordinate newtypes for old and new assemblies.
//!
//! All intervals are 0-based half-open. The types are deliberately distinct so
//! that mixing old and new coordinates is a compile error. They stay minimal:
//! constructors, `len`, `contains` / `overlaps`. No generic trait machinery.

/// A position on the old assembly.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct OldPos(pub u64);

/// A position on the new assembly.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct NewPos(pub u64);

impl OldPos {
    #[inline]
    pub fn get(self) -> u64 {
        self.0
    }
}

impl NewPos {
    #[inline]
    pub fn get(self) -> u64 {
        self.0
    }
}

/// A half-open interval `[start, end)` in old coordinates.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct OldIv {
    pub start: OldPos,
    pub end: OldPos,
}

/// A half-open interval `[start, end)` in new coordinates.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct NewIv {
    pub start: NewPos,
    pub end: NewPos,
}

impl OldIv {
    #[inline]
    pub fn new(start: u64, end: u64) -> Self {
        debug_assert!(start <= end);
        Self { start: OldPos(start), end: OldPos(end) }
    }

    #[inline]
    pub fn len(&self) -> u64 {
        self.end.0 - self.start.0
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.start.0 == self.end.0
    }

    /// True if `p` lies inside the interval (half-open).
    #[inline]
    pub fn contains(&self, p: OldPos) -> bool {
        self.start.0 <= p.0 && p.0 < self.end.0
    }

    /// True if `other` is fully contained in `self`.
    #[inline]
    pub fn contains_iv(&self, other: OldIv) -> bool {
        self.start.0 <= other.start.0 && other.end.0 <= self.end.0
    }

    /// True if the two intervals share at least one base.
    #[inline]
    pub fn overlaps(&self, other: OldIv) -> bool {
        self.start.0 < other.end.0 && other.start.0 < self.end.0
    }
}

impl NewIv {
    #[inline]
    pub fn new(start: u64, end: u64) -> Self {
        debug_assert!(start <= end);
        Self { start: NewPos(start), end: NewPos(end) }
    }

    #[inline]
    pub fn len(&self) -> u64 {
        self.end.0 - self.start.0
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.start.0 == self.end.0
    }

    #[inline]
    pub fn contains(&self, p: NewPos) -> bool {
        self.start.0 <= p.0 && p.0 < self.end.0
    }

    #[inline]
    pub fn overlaps(&self, other: NewIv) -> bool {
        self.start.0 < other.end.0 && other.start.0 < self.end.0
    }
}
