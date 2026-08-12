# Interactive thesis demo

A responsive browser reconstruction of the shrinking-ball medial-axis algorithm used in the native BEng thesis project.

## Development

```bash
npm ci
npm run dev
```

Open the local URL printed by the development server. The page is self-contained and does not require a database or external API.

## Quality checks

```bash
npm run lint
npm test
```

`npm test` creates a production build and verifies the server-rendered title, project copy, accessible canvas description and interactive controls.

## Notes

- The browser visualisation is a two-dimensional educational reconstruction, not a WebAssembly build of PCL.
- The exact native implementation and NumPy workflow live one directory above this app.
- Social metadata resolves its absolute image URL from the incoming request host.
