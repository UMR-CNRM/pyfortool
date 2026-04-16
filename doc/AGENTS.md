# AGENTS.md – PyForTool Documentation

This file provides specific guidance for working with PyForTool's documentation.

## Scope

This file applies to all files in the `doc/` directory:
- `Documentation.md` - User guide
- `doc/developer/` - Developer documentation
- `doc/doxygen/` - Doxygen configuration and mainpage

## Key Rule: Always Use .md Extensions

**Doxygen automatically transforms .md links to .html.** Using .html extensions in documentation source files will result in broken links.

| ✅ Correct | ❌ Incorrect |
|-----------|--------------|
| `doc/developer/architecture.md` | `doc/developer/architecture.html` |
| `doc/Documentation.md` | `doc/Documentation.html` |
| `concepts.md` | `concepts.html` |

## Link Examples

### Absolute paths from project root
```markdown
[Architecture Guide](../doc/developer/architecture.md)
[User Guide](../doc/Documentation.md)
```

### Relative paths within doc/
```markdown
[Architecture Guide](developer/architecture.md)
[Core Concepts](developer/concepts.md)
[Testing Guide](developer/testing.md)
```

### From Doxygen mainpage (doc/doxygen/mainpage.md)
```markdown
[User's Guide](../Documentation.md)
[Architecture Guide](../developer/architecture.md)
```

## Common Mistakes to Avoid

1. **Using .html extensions:**
   ```markdown
   ❌ [Architecture Guide](md__home_sriette_GIT_pyfortool_doc_developer_architecture.html)
   ✅ [Architecture Guide](developer/architecture.md)
   ```

2. **Hardcoded paths with usernames:**
   ```markdown
   ❌ md__home_sriette_GIT_pyfortool_doc_developer_architecture.html
   ✅ doc/developer/architecture.md
   ```

3. **Forgetting to update links when moving files:**
   - Always verify links work after reorganizing docs
   - Run Doxygen and check generated HTML for broken links

## Testing Links

After modifying documentation:
1. Build Doxygen: `cd doc/doxygen && doxygen doxygen_config`
2. Check generated HTML files in `doc/html/`
3. Verify cross-references work in browser

## Related Files

- Root AGENTS.md: `../AGENTS.md` - General project guidance
- User guide: `Documentation.md` - End-user documentation
- Developer docs: `developer/` - Architecture, concepts, testing