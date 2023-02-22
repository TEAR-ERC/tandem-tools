--
-- Print the contents of `tbl` (recursively) with the indentation depth set by `indent`.
-- Example usage:
--
--   require("table_print")  
--
--   volume_properties = {}
--   volume_properties["crust"] = { cs = 3400.0, rho = 2700.0 }
--   volume_properties["mantle"] = { cs = 3000.0, rho = 3200.0 }
--   tprint(volume_properties)
-- Contact: dmay@ucsd.edu
--
function tprint(tbl, indent)
  if not indent then indent = 0 end

  for k, v in pairs(tbl) do
    formatting  = string.rep(" ", indent) .. tostring(k) .. ": "
    if type(v) == "table" then
      print(formatting)
      tprint(v, indent+1)
    else
      print(formatting .. tostring(v))
    end
  end
end
