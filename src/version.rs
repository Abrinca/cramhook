//! Output format version.
//!
//! Every machine-readable output written by the tool carries the same
//! `format_version` so consumers can detect breaking changes. It is the first
//! key of every JSON document (`format_version` in JSON) and a VCF header meta
//! line (`##format_version=1`).
pub const FORMAT_VERSION: u32 = 1;

/// Wraps a serializable output document so that `format_version` is written
/// as the first JSON key, followed by the document's own fields.
#[derive(serde::Serialize)]
pub struct Versioned<'a, T: serde::Serialize> {
    pub format_version: u32,
    #[serde(flatten)]
    pub value: &'a T,
}

/// Wrap `value` for writing with the current [`FORMAT_VERSION`] first.
pub fn versioned<T: serde::Serialize>(value: &T) -> Versioned<'_, T> {
    Versioned { format_version: FORMAT_VERSION, value }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_version_is_the_first_key() {
        let doc = serde_json::json!({"b": 1, "a": 2});
        let s = serde_json::to_string(&versioned(&doc)).unwrap();
        assert_eq!(s, "{\"format_version\":1,\"a\":2,\"b\":1}");
    }
}
