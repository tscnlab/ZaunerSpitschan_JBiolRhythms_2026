-- Word-only presentation of the existing manuscript metadata and section levels.
-- Run after the main MeLiDos paper's authors-block filter.
function Pandoc(doc)
  if not FORMAT:match('docx') then return doc end

  doc = doc:walk({Link = function(link)
    -- Word travels separately from the website. Named supplements point to the
    -- corresponding published supplement, including section anchors.
    if link.target:match('%.qmd') then
      link.target = 'https://tscnlab.github.io/ZaunerSpitschan_JBiolRhythms_2026/' .. link.target:gsub('%.qmd', '.html')
    end
    return link
  end, Header = function(h)
    h.level = math.max(1, h.level - 1)
    return h
  end})

  local abstract = doc.meta.abstract
  local front = pandoc.List()
  if abstract then
    front:insert(pandoc.Header(1, 'Abstract'))
    front:extend(abstract)
    doc.meta.abstract = nil
  end
  if doc.meta.keywords then
    local words = {}
    for _, keyword in ipairs(doc.meta.keywords) do
      table.insert(words, pandoc.utils.stringify(keyword))
    end
    local inlines = pandoc.List({pandoc.Strong('Keywords:'), pandoc.Space()})
    inlines:extend(pandoc.Inlines(table.concat(words, '; ')))
    front:insert(pandoc.Para(inlines))
  end

  local body = pandoc.List()
  local inserted = false
  for _, block in ipairs(doc.blocks) do
    if not inserted and block.t == 'Header' then
      body:extend(front)
      inserted = true
    end
    body:insert(block)
  end
  if not inserted then body:extend(front) end
  doc.blocks = body
  return doc
end
