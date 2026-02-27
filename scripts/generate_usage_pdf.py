from pathlib import Path
import textwrap

def escape_pdf_text(s: str) -> str:
    return s.replace('\\', '\\\\').replace('(', '\\(').replace(')', '\\)')

def build_page_content(lines, page_w=595, page_h=842, margin=50, font_size=10, leading=14):
    y = page_h - margin
    content = ['BT', f'/F1 {font_size} Tf']
    for line in lines:
        if y < margin:
            break
        content.append(f'1 0 0 1 {margin} {y} Tm ({escape_pdf_text(line)}) Tj')
        y -= leading
    content.append('ET')
    return '\n'.join(content)

def write_pdf(pages, out_path: Path):
    objs = []

    # 1: Catalog, 2: Pages, 3: Font
    objs.append('<< /Type /Catalog /Pages 2 0 R >>')
    objs.append('<< /Type /Pages /Kids [ ' + ' '.join(f'{4+i*2} 0 R' for i in range(len(pages))) + f' ] /Count {len(pages)} >>')
    objs.append('<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>')

    for i, content in enumerate(pages):
        page_obj_num = 4 + i * 2
        content_obj_num = page_obj_num + 1
        objs.append(f'<< /Type /Page /Parent 2 0 R /MediaBox [0 0 595 842] /Resources << /Font << /F1 3 0 R >> >> /Contents {content_obj_num} 0 R >>')
        stream = content.encode('latin-1', errors='ignore')
        objs.append(f'<< /Length {len(stream)} >>\nstream\n' + stream.decode('latin-1', errors='ignore') + '\nendstream')

    pdf = ['%PDF-1.4']
    offsets = [0]
    cur = len(pdf[0]) + 1

    for idx, obj in enumerate(objs, start=1):
        offsets.append(cur)
        chunk = f'{idx} 0 obj\n{obj}\nendobj\n'
        pdf.append(chunk)
        cur += len(chunk)

    xref_pos = cur
    xref = [f'xref\n0 {len(objs)+1}', '0000000000 65535 f ']
    for off in offsets[1:]:
        xref.append(f'{off:010d} 00000 n ')
    trailer = f'trailer\n<< /Size {len(objs)+1} /Root 1 0 R >>\nstartxref\n{xref_pos}\n%%EOF\n'

    out = '\n'.join(pdf) + '\n' + '\n'.join(xref) + '\n' + trailer
    out_path.write_bytes(out.encode('latin-1', errors='ignore'))


def main():
    src = Path('docs/usage_guide.md')
    dst = Path('docs/usage_guide.pdf')

    raw = src.read_text(encoding='utf-8')
    # PDF core writer uses Helvetica latin-1; strip non-latin chars to keep PDF valid.
    ascii_text = raw.encode('latin-1', errors='ignore').decode('latin-1')
    lines = []
    for line in ascii_text.splitlines():
        if line.strip().startswith('```'):
            lines.append('--- code ---')
            continue
        wrapped = textwrap.wrap(line, width=95) if line else ['']
        lines.extend(wrapped)

    page_capacity = 52
    pages = []
    for i in range(0, len(lines), page_capacity):
        pages.append(build_page_content(lines[i:i+page_capacity]))

    write_pdf(pages, dst)
    print(f'Wrote {dst} ({len(pages)} pages)')

if __name__ == '__main__':
    main()
